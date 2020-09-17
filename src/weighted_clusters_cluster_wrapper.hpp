#pragma once

#include "weighted_clusters_iterator.hpp"

#include <type_traits>

class WeightedClusters;
template <typename>
struct WeightedClustersCluster;

template <typename T>
class WeightedClustersClusterWrapper {
  using weighted_clusters_type =
      std::conditional_t<std::is_const<T>::value, const WeightedClusters,
                         WeightedClusters>;
  template <typename>
  friend class WeightedClustersClusterIterator;

public:
  using concrete_type = WeightedClustersCluster<std::decay_t<T>>;
  using iterator = WeightedClustersIterator<T, false>;
  using reference = std::conditional_t<std::is_rvalue_reference_v<T>, T, T&>;

  WeightedClustersClusterWrapper() = default;
  WeightedClustersClusterWrapper(weighted_clusters_type& weightedClusters,
                                 std::ptrdiff_t clusterIndex) noexcept;

  WeightedClustersClusterWrapper& operator=(concrete_type const& cluster) const
      noexcept(std::is_nothrow_copy_assignable_v<std::decay_t<T>>);
  WeightedClustersClusterWrapper& operator=(concrete_type&& cluster) const
      noexcept(std::is_nothrow_move_assignable_v<std::decay_t<T>>);

  iterator begin() const noexcept;
  iterator end() const noexcept;

  reference operator[](std::size_t index) const noexcept;

  std::size_t size() const noexcept;

  operator concrete_type() const noexcept(false);
  std::size_t index() const noexcept;

private:
  weighted_clusters_type* weightedClusters;
  std::ptrdiff_t clusterIndex;
};

#include "weighted_clusters_cluster_wrapper_impl.hpp"

#pragma once

#include "weighted_clusters_cluster_wrapper.hpp"

#include <range/v3/core.hpp>
#include <type_traits>

class WeightedClusters;
template <typename> struct WeightedClustersCluster;

template <typename T> class WeightedClustersClusterIterator {
  using weighted_clusters_type =
      std::conditional_t<std::is_const<T>::value, const WeightedClusters,
                         WeightedClusters>;

public:
  using difference_type = std::ptrdiff_t;
  using wrapper_type = WeightedClustersClusterWrapper<T>;
  using value_type = WeightedClustersCluster<std::decay_t<T>>;
  using reference = wrapper_type;
  using iterator_category = ranges::random_access_iterator_tag;

  WeightedClustersClusterIterator() = default;
  explicit WeightedClustersClusterIterator(
      weighted_clusters_type &weightedClusters) noexcept;
  WeightedClustersClusterIterator(weighted_clusters_type &weightedClusters,
                                  difference_type clusterIndex) noexcept;

  WeightedClustersClusterIterator &operator++() noexcept;
  WeightedClustersClusterIterator operator++(int) noexcept;
  WeightedClustersClusterIterator &operator--() noexcept;
  WeightedClustersClusterIterator operator--(int) noexcept;
  WeightedClustersClusterIterator &operator+=(difference_type offset) noexcept;
  WeightedClustersClusterIterator &operator-=(difference_type offset) noexcept;
  WeightedClustersClusterIterator
  operator+(difference_type offset) const noexcept;
  WeightedClustersClusterIterator
  operator-(difference_type offset) const noexcept;

  template <typename U>
  friend WeightedClustersClusterIterator<U>
  operator+(typename WeightedClustersClusterIterator<U>::difference_type offset,
            const WeightedClustersClusterIterator<U> &iter) noexcept;

  reference operator*() const noexcept;
  reference operator[](difference_type offset) const noexcept;

  bool operator==(const WeightedClustersClusterIterator &other) const noexcept;
  bool operator!=(const WeightedClustersClusterIterator &other) const noexcept;
  bool operator<(const WeightedClustersClusterIterator &other) const noexcept;
  bool operator<=(const WeightedClustersClusterIterator &other) const noexcept;
  bool operator>=(const WeightedClustersClusterIterator &other) const noexcept;
  bool operator>(const WeightedClustersClusterIterator &other) const noexcept;
  difference_type
  operator-(const WeightedClustersClusterIterator &other) const noexcept;

private:
  wrapper_type wrapper;
};

#include "weighted_clusters_cluster_iterator_impl.hpp"

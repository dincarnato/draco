#pragma once

#include "clusters_traits.hpp"

#include <type_traits>

template <typename, bool>
struct HardClusterWrapper;
template <typename, bool>
struct HardClustersIterator;

template <typename Cluster, bool complemented>
struct HardClustersWrapper {
  using cluster_type = typename HardClustersTraits<Cluster>::cluster_type;
  using index_type = typename HardClustersTraits<Cluster>::index_type;
  using cluster_wrapper = HardClusterWrapper<Cluster, complemented>;
  using iterator = HardClustersIterator<Cluster, complemented>;

  HardClustersWrapper() = default;
  explicit HardClustersWrapper(Cluster& clusters) noexcept;

  cluster_wrapper operator[](index_type index) const noexcept;

  iterator begin() const noexcept;
  iterator end() const noexcept;

  std::size_t size() const noexcept;
  std::size_t elements_size() const noexcept;

private:
  Cluster* _clusters;
};

#include "hard_clusters_wrapper_impl.hpp"

#pragma once

#include "clusters_traits.hpp"
#include <type_traits>

template <typename Cluster, bool complemented>
struct HardClustersElementWrapper {
  using cluster_type = typename HardClustersTraits<Cluster>::cluster_type;
  using index_type = typename HardClustersTraits<Cluster>::index_type;
  using element_type = std::conditional_t<std::is_const_v<Cluster>,
                                          const index_type, index_type>;

  HardClustersElementWrapper() = default;
  HardClustersElementWrapper(Cluster &cluster, std::size_t index) noexcept;

  index_type get() const noexcept;
  void set(index_type clusterIndex) const noexcept;
  void clear() const noexcept;
  bool isCluster(index_type index) const noexcept;

private:
  element_type *_element;
#ifndef NDEBUG
  Cluster *_cluster;
#endif
};

#include "hard_clusters_element_wrapper_impl.hpp"

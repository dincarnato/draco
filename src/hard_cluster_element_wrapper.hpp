#pragma once

#include "clusters_traits.hpp"
#include <type_traits>

template <typename, bool>
struct HardClusterWrapper;

template <typename Cluster, bool complemented>
struct HardClusterElementWrapper {
  using cluster_type = typename HardClustersTraits<Cluster>::cluster_type;
  using index_type = typename HardClustersTraits<Cluster>::index_type;
  using element_type = std::conditional_t<std::is_const_v<Cluster>,
                                          const index_type, index_type>;
  static constexpr bool isComplemented = complemented;

  template <typename, bool>
  friend struct HardClusterElementIterator;

  HardClusterElementWrapper() = default;
  template <typename OtherCluster>
  HardClusterElementWrapper(
      const HardClusterWrapper<OtherCluster, complemented>& clusterWrapper,
      std::size_t elementIndex) noexcept;

  HardClusterElementWrapper const& operator=(bool value) const noexcept;

  bool get() const noexcept;
  void set() const noexcept;
  void clear() const noexcept;
  operator bool() const noexcept;
  void swap(HardClusterElementWrapper other) const noexcept;
  void swapMatched(HardClusterElementWrapper other) const noexcept;

private:
  element_type* _element;
  index_type clusterIndex;
#ifndef NDEBUG
  Cluster* _cluster;
#endif
};

#include "hard_cluster_element_wrapper_impl.hpp"

#pragma once

#include "clusters_traits.hpp"

#include <type_traits>

struct HardCluster;
template <typename, bool> struct HardClusterElementWrapper;
template <typename, bool> struct HardClustersIterator;
template <typename, bool> struct HardClusterElementIterator;

template <typename Cluster, bool complemented> struct HardClusterWrapper {
  template <typename, bool> friend struct HardClusterElementWrapper;
  template <typename, bool> friend struct HardClustersIterator;

  using concrete_type = HardCluster;
  using cluster_type = typename HardClustersTraits<Cluster>::cluster_type;
  using index_type = typename HardClustersTraits<Cluster>::index_type;
  using element_wrapper = HardClusterElementWrapper<Cluster, complemented>;
  using iterator = HardClusterElementIterator<Cluster, complemented>;
  static constexpr bool isComplemented = complemented;

  HardClusterWrapper() = default;
  HardClusterWrapper(Cluster &cluster, index_type index) noexcept;

  explicit operator concrete_type() const noexcept(false);
  HardClusterWrapper const &operator=(concrete_type const &rhs) const noexcept;
  HardClusterWrapper const &operator=(concrete_type &&rhs) const noexcept;

  element_wrapper operator[](std::size_t index) const noexcept;
  void set(std::size_t index) const noexcept;
  void clear(std::size_t index) const noexcept;
  void moveTo(std::size_t index, std::size_t otherClusterIndex) const noexcept;

  iterator begin() const noexcept;
  iterator end() const noexcept;

  index_type index() const noexcept;
  std::size_t size() const noexcept;

private:
  Cluster *_cluster;
  index_type _index;
};

#include "hard_cluster_wrapper_impl.hpp"

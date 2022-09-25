#pragma once

#include "clusters_traits.hpp"

template <typename, bool> struct HardClustersWrapper;
template <typename, bool> struct HardClusterWrapper;
template <typename, bool> struct HardClustersElementWrapper;

template <typename Base> struct HardClustersBaseComplement {
  using hard_clusters_base_type = Base;
  using index_type = typename HardClustersTraits<Base>::index_type;
  using clusters_wrapper_type =
      HardClustersWrapper<hard_clusters_base_type, true>;
  using cluster_wrapper_type =
      HardClusterWrapper<hard_clusters_base_type, true>;
  using element_wrapper =
      HardClustersElementWrapper<hard_clusters_base_type, true>;

  template <typename, bool> friend struct HardClusterWrapper;

  explicit HardClustersBaseComplement(
      hard_clusters_base_type &clusters) noexcept;

  cluster_wrapper_type cluster(index_type index) const noexcept;
  clusters_wrapper_type clusters() const noexcept;
  element_wrapper operator[](std::size_t index) const noexcept;
  hard_clusters_base_type &complement() const noexcept;

private:
  hard_clusters_base_type *_clusters;
};

#include "hard_clusters_base_complement_impl.hpp"

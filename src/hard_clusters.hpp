#pragma once

#include "clusters_traits.hpp"
#include "hard_cluster.hpp"
#include "hard_cluster_wrapper.hpp"
#include "hard_clusters_base_complement.hpp"
#include "hard_clusters_element_wrapper.hpp"
#include "hard_clusters_iterator.hpp"
#include "hard_clusters_wrapper.hpp"

#include <range/v3/iterator/concepts.hpp>

#include <vector>

template <typename Type>
struct HardClustersBase {
  using index_type = Type;
  using clusters_type = std::vector<Type>;
  using clusters_wrapper_type = HardClustersWrapper<HardClustersBase, false>;
  using const_clusters_wrapper_type =
      HardClustersWrapper<const HardClustersBase, false>;
  using cluster_wrapper_type = HardClusterWrapper<HardClustersBase, false>;
  using const_cluster_wrapper_type =
      HardClusterWrapper<const HardClustersBase, false>;
  using element_wrapper = HardClustersElementWrapper<HardClustersBase, false>;
  using const_element_wrapper =
      HardClustersElementWrapper<const HardClustersBase, false>;
  using complement_type = HardClustersBaseComplement<HardClustersBase>;
  using const_complement_type =
      HardClustersBaseComplement<const HardClustersBase>;
  using iterator = typename clusters_type::iterator;
  using const_iterator = typename clusters_type::const_iterator;

  template <typename, bool>
  friend struct HardClusterWrapper;
  template <typename, bool>
  friend struct HardClustersWrapper;
  template <typename, bool>
  friend struct HardClusterElementWrapper;
  template <typename, bool>
  friend struct HardClustersElementWrapper;
  template <typename>
  friend struct HardClustersBaseComplement;

  HardClustersBase() = default;
  HardClustersBase(std::size_t elements, Type clusters) noexcept;
  HardClustersBase(std::size_t elements, Type clusters,
                   Type defaultCluster) noexcept;

  cluster_wrapper_type cluster(index_type index);
  const_cluster_wrapper_type cluster(index_type index) const;

  clusters_wrapper_type clusters();
  const_clusters_wrapper_type clusters() const;

  element_wrapper operator[](std::size_t index);
  const_element_wrapper operator[](std::size_t index) const;

  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;

  void clear();

  complement_type complement();
  const_complement_type complement() const;

  std::size_t getElementsSize() const;
  index_type getClustersSize() const;

  void extendMinorCluster(std::size_t maxDistance = 10);
  void removeSmallSpans(std::size_t minSize);

private:
  Type nClusters = 0;
  clusters_type _clusters;
};

#include "hard_clusters_impl.hpp"

using HardClusters = HardClustersBase<unsigned char>;
static_assert(is_clusters_v<HardClusters>);

static_assert(ranges::RandomAccessIterator<typename HardClusters::iterator>);
static_assert(
    ranges::RandomAccessIterator<typename HardClusters::const_iterator>);

static_assert(ranges::RandomAccessIterator<
              typename HardClusters::cluster_wrapper_type::iterator>);
static_assert(ranges::RandomAccessIterator<
              typename HardClusters::const_cluster_wrapper_type::iterator>);

static_assert(ranges::RandomAccessIterator<
              typename HardClusters::clusters_wrapper_type::iterator>);
static_assert(ranges::RandomAccessIterator<
              typename HardClusters::const_clusters_wrapper_type::iterator>);

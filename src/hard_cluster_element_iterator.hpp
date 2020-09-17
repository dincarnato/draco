#pragma once

#include "clusters_traits.hpp"
#include <range/v3/core.hpp>
#include <type_traits>

template <typename, bool>
struct HardClusterWrapper;
template <typename, bool>
struct HardClusterElementWrapper;

template <typename Cluster, bool complemented>
struct HardClusterElementIterator {
  using cluster_type = typename HardClustersTraits<Cluster>::cluster_type;
  using index_type = typename HardClustersTraits<Cluster>::index_type;
  using hard_cluster_wrapper_type = HardClusterWrapper<Cluster, complemented>;
  using wrapper_type = HardClusterElementWrapper<Cluster, complemented>;
  using value_type = bool;
  using reference = wrapper_type;
  using difference_type = std::ptrdiff_t;
  using iterator_category = ranges::random_access_iterator_tag;

  HardClusterElementIterator() = default;
  HardClusterElementIterator(hard_cluster_wrapper_type const& cluster,
                             std::size_t elementIndex) noexcept;

private:
  using self = HardClusterElementIterator;

public:
  self& operator++() noexcept;
  self operator++(int) noexcept;
  self& operator--() noexcept;
  self operator--(int) noexcept;
  self& operator+=(difference_type offset) noexcept;
  self& operator-=(difference_type offset) noexcept;
  self operator+(difference_type offset) const noexcept;
  self operator-(difference_type offset) const noexcept;

  template <typename _Cluster, bool _complemented>
  friend HardClusterElementIterator<_Cluster, _complemented> operator+(
      typename HardClusterElementIterator<
          _Cluster, _complemented>::difference_type offset,
      const HardClusterElementIterator<_Cluster, _complemented>& iter) noexcept;

  bool operator<(const self& other) const noexcept;
  bool operator>(const self& other) const noexcept;
  bool operator>=(const self& other) const noexcept;
  bool operator<=(const self& other) const noexcept;
  bool operator==(const self& other) const noexcept;
  bool operator!=(const self& other) const noexcept;
  difference_type operator-(const self& other) const noexcept;

  reference operator*() const noexcept;
  reference operator[](difference_type offset) const noexcept;

private:
  wrapper_type wrapper;
};

#include "hard_cluster_element_iterator_impl.hpp"

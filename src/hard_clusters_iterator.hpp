#pragma once

#include "hard_cluster_wrapper.hpp"

#include <cstdint>
#include <range/v3/core.hpp>
#include <type_traits>

struct HardCluster;

template <typename Cluster, bool complemented>
struct HardClustersIterator {
  using cluster_type = typename HardClustersTraits<Cluster>::cluster_type;
  using index_type = typename HardClustersTraits<Cluster>::index_type;
  using iterator_category = ranges::random_access_iterator_tag;
  using wrapper_type = HardClusterWrapper<Cluster, complemented>;
  using value_type = HardCluster;
  using reference = wrapper_type;
  using difference_type = std::conditional_t<
      sizeof(index_type) == 1, std::int16_t,
      std::conditional_t<sizeof(index_type) == 2, std::int32_t,
                         std::conditional_t<sizeof(index_type) == 4,
                                            std::int64_t, std::ptrdiff_t>>>;

  HardClustersIterator() = default;
  explicit HardClustersIterator(Cluster& cluster,
                                index_type offset = 0) noexcept;

  HardClustersIterator const& operator=(value_type const& rhs) const noexcept;
  HardClustersIterator const& operator=(value_type&& rhs) const noexcept;

  bool operator<(const HardClustersIterator& other) const noexcept;
  bool operator<=(const HardClustersIterator& other) const noexcept;
  bool operator>=(const HardClustersIterator& other) const noexcept;
  bool operator>(const HardClustersIterator& other) const noexcept;
  bool operator==(const HardClustersIterator& other) const noexcept;
  bool operator!=(const HardClustersIterator& other) const noexcept;
  difference_type operator-(const HardClustersIterator& other) const noexcept;

  HardClustersIterator& operator++() noexcept;
  HardClustersIterator operator++(int) noexcept;
  HardClustersIterator& operator--() noexcept;
  HardClustersIterator operator--(int) noexcept;
  HardClustersIterator& operator+=(difference_type offset) noexcept;
  HardClustersIterator& operator-=(difference_type offset) noexcept;
  HardClustersIterator operator+(difference_type offset) const noexcept;
  HardClustersIterator operator-(difference_type offset) const noexcept;

  template <typename _Cluster, bool _complemented>
  friend HardClustersIterator<_Cluster, _complemented> operator+(
      typename HardClustersIterator<_Cluster, _complemented>::difference_type,
      const HardClustersIterator<_Cluster, _complemented>& iter) noexcept;

  reference operator*() const noexcept;
  reference operator[](difference_type offset) const noexcept;

private:
  wrapper_type wrapper;
};

#include "hard_clusters_iterator_impl.hpp"

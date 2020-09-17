#pragma once

#include <iterator>
#include <range/v3/core.hpp>
#include <type_traits>

class WeightedClusters;
template <typename>
struct WeightedClustersConcreteSpan;
template <typename>
class WeightedClustersSpan;

template <typename T, bool useSpans = true>
class WeightedClustersIterator {
  using weighted_clusters_type =
      std::conditional_t<std::is_const<T>::value, const WeightedClusters,
                         WeightedClusters>;
  using weighted_clusters_span_type = WeightedClustersSpan<T>;

  struct EndTag {};

public:
  using difference_type = std::ptrdiff_t;
  using value_type = std::conditional_t<
      useSpans, WeightedClustersConcreteSpan<std::decay_t<T>>, std::decay_t<T>>;
  using reference = std::conditional_t<
      useSpans, weighted_clusters_span_type,
      std::conditional_t<std::is_rvalue_reference_v<T>, T, T&>>;
  using iterator_category = ranges::random_access_iterator_tag;

  static constexpr EndTag endTag{};

  explicit WeightedClustersIterator() = default;
  explicit WeightedClustersIterator(weighted_clusters_type& weightedClusters,
                                    difference_type offset = 0) noexcept;
  WeightedClustersIterator(weighted_clusters_type& weightedClusters,
                           difference_type offset, EndTag) noexcept;

  template <bool _useSpans = useSpans,
            typename = std::enable_if_t<not _useSpans>>
  WeightedClustersIterator(weighted_clusters_type& weightedClusters,
                           difference_type clusterOffset,
                           difference_type elementOffset) noexcept;

  WeightedClustersIterator& operator++() noexcept;
  WeightedClustersIterator operator++(int) noexcept;
  WeightedClustersIterator& operator--() noexcept;
  WeightedClustersIterator operator--(int) noexcept;
  WeightedClustersIterator& operator+=(difference_type offset) noexcept;
  WeightedClustersIterator& operator-=(difference_type offset) noexcept;
  WeightedClustersIterator operator+(difference_type offset) const noexcept;
  WeightedClustersIterator operator-(difference_type offset) const noexcept;

  template <typename U, bool _useSpans>
  friend WeightedClustersIterator<U, _useSpans> operator+(
      typename WeightedClustersIterator<U, _useSpans>::difference_type offset,
      const WeightedClustersIterator<U, _useSpans>& iter) noexcept;

  reference operator*() const noexcept;
  reference operator[](difference_type index) const noexcept;

  bool operator<(const WeightedClustersIterator& other) const noexcept;
  bool operator>(const WeightedClustersIterator& other) const noexcept;
  bool operator>=(const WeightedClustersIterator& other) const noexcept;
  bool operator<=(const WeightedClustersIterator& other) const noexcept;
  bool operator==(const WeightedClustersIterator& other) const noexcept;
  bool operator!=(const WeightedClustersIterator& other) const noexcept;
  difference_type operator-(const WeightedClustersIterator& other) const
      noexcept;

private:
  weighted_clusters_span_type _data;
  std::size_t clusters;
};

#include "weighted_clusters_iterator_impl.hpp"

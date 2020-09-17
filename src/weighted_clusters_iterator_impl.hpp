#pragma once

#include "weighted_clusters.hpp"
#include "weighted_clusters_concrete_span.hpp"
#include "weighted_clusters_iterator.hpp"
#include "weighted_clusters_span.hpp"

#include <cassert>

template <typename T, bool useSpans>
WeightedClustersIterator<T, useSpans>::WeightedClustersIterator(
    weighted_clusters_type& weightedClusters, difference_type offset) noexcept
    : _data([&] {
        assert(offset >= 0);

        if constexpr (useSpans) {
          // offset is an element offset
          assert(offset <=
                 static_cast<difference_type>(weightedClusters.elements));
          return weighted_clusters_span_type(
              weightedClusters.weights, weightedClusters._clusters,
              offset *
                  static_cast<difference_type>(weightedClusters._clusters));
        } else {
          // offset is a cluster offset
          assert(offset <=
                 static_cast<difference_type>(weightedClusters._clusters));
          return weighted_clusters_span_type(weightedClusters.weights,
                                             weightedClusters.elements, offset);
        }
      }()),
      clusters(weightedClusters._clusters) {}

template <typename T, bool useSpans>
WeightedClustersIterator<T, useSpans>::WeightedClustersIterator(
    weighted_clusters_type& weightedClusters, difference_type offset,
    EndTag) noexcept
    : _data([&] {
        assert(offset >= 0);

        if constexpr (useSpans) {
          // offset is an element offset
          assert(offset <= weightedClusters._elements);
          return weighted_clusters_span_type(
              weightedClusters.weights, weightedClusters._clusters,
              (offset + 1) *
                  static_cast<difference_type>(weightedClusters._clusters));
        } else {
          // offset is a cluster offset
          assert(offset <=
                 static_cast<difference_type>(weightedClusters._clusters));
          return weighted_clusters_span_type(
              weightedClusters.weights, weightedClusters.elements,
              offset + static_cast<difference_type>(
                           weightedClusters.weights.size()));
        }
      }()),
      clusters(weightedClusters._clusters) {}

template <typename T, bool useSpans>
template <bool, typename>
WeightedClustersIterator<T, useSpans>::WeightedClustersIterator(
    weighted_clusters_type& weightedClusters, difference_type clusterOffset,
    difference_type elementOffset) noexcept
    : _data(weightedClusters.weights, weightedClusters.elements,
            elementOffset *
                    static_cast<difference_type>(weightedClusters._clusters) +
                clusterOffset),
      clusters(weightedClusters._clusters) {
  assert(clusterOffset >= 0);
  assert(clusterOffset <=
         static_cast<difference_type>(weightedClusters._clusters));
  assert(elementOffset >= 0);
  assert(elementOffset <=
         static_cast<difference_type>(weightedClusters.elements));
}

template <typename T, bool useSpans>
auto
WeightedClustersIterator<T, useSpans>::operator++() noexcept
    -> WeightedClustersIterator& {
  _data.index += static_cast<difference_type>(clusters);
  return *this;
}

template <typename T, bool useSpans>
auto
WeightedClustersIterator<T, useSpans>::operator++(int) noexcept
    -> WeightedClustersIterator {
  auto other = *this;
  _data.index += static_cast<difference_type>(clusters);
  return other;
}

template <typename T, bool useSpans>
auto
WeightedClustersIterator<T, useSpans>::operator--() noexcept
    -> WeightedClustersIterator& {
  _data.index -= static_cast<difference_type>(clusters);
  return *this;
}

template <typename T, bool useSpans>
auto
WeightedClustersIterator<T, useSpans>::operator--(int) noexcept
    -> WeightedClustersIterator {
  auto other = *this;
  _data.index -= static_cast<difference_type>(clusters);
  return other;
}

template <typename T, bool useSpans>
auto
WeightedClustersIterator<T, useSpans>::
operator+=(difference_type offset) noexcept -> WeightedClustersIterator& {
  _data.index += static_cast<difference_type>(clusters) * offset;
  return *this;
}

template <typename T, bool useSpans>
auto
WeightedClustersIterator<T, useSpans>::
operator-=(difference_type offset) noexcept -> WeightedClustersIterator& {
  _data.index -= static_cast<difference_type>(clusters) * offset;
  return *this;
}

template <typename T, bool useSpans>
auto
WeightedClustersIterator<T, useSpans>::operator+(difference_type offset) const
    noexcept -> WeightedClustersIterator {
  auto other = *this;
  other._data.index += static_cast<difference_type>(clusters) * offset;
  return other;
}

template <typename T, bool useSpans>
auto
WeightedClustersIterator<T, useSpans>::operator-(difference_type offset) const
    noexcept -> WeightedClustersIterator {
  auto other = *this;
  other._data.index -= static_cast<difference_type>(clusters) * offset;

  return other;
}

template <typename T, bool useSpans>
auto WeightedClustersIterator<T, useSpans>::operator*() const noexcept
    -> reference {
  if constexpr (useSpans)
    return _data;
  else
    return *_data;
}

template <typename T, bool useSpans>
auto WeightedClustersIterator<T, useSpans>::
operator[](difference_type index) const noexcept -> reference {
  if constexpr (useSpans) {
    auto newData = _data;
    newData.index += static_cast<difference_type>(clusters) * index;
    return newData;
  } else
    return _data[static_cast<std::size_t>(
        index * static_cast<difference_type>(clusters))];
}

template <typename T, bool useSpans>
bool
WeightedClustersIterator<T, useSpans>::
operator<(const WeightedClustersIterator& other) const noexcept {
  return _data.index < other._data.index;
}

template <typename T, bool useSpans>
bool
WeightedClustersIterator<T, useSpans>::
operator>(const WeightedClustersIterator& other) const noexcept {
  return _data.index > other._data.index;
}

template <typename T, bool useSpans>
bool
WeightedClustersIterator<T, useSpans>::
operator>=(const WeightedClustersIterator& other) const noexcept {
  return _data.index >= other._data.index;
}

template <typename T, bool useSpans>
bool
WeightedClustersIterator<T, useSpans>::
operator<=(const WeightedClustersIterator& other) const noexcept {
  return _data._data <= other._data._data;
}

template <typename T, bool useSpans>
bool
WeightedClustersIterator<T, useSpans>::
operator==(const WeightedClustersIterator& other) const noexcept {
  return _data.index == other._data.index;
}

template <typename T, bool useSpans>
bool
WeightedClustersIterator<T, useSpans>::
operator!=(const WeightedClustersIterator& other) const noexcept {
  return _data.index != other._data.index;
}

template <typename T, bool useSpans>
auto
WeightedClustersIterator<T, useSpans>::
operator-(const WeightedClustersIterator& other) const noexcept
    -> difference_type {
  assert(clusters == other.clusters);
  return static_cast<difference_type>(_data.index - other._data.index) /
         static_cast<difference_type>(clusters);
}

template <typename T, bool useSpans>
WeightedClustersIterator<T, useSpans>
operator+(
    typename WeightedClustersIterator<T, useSpans>::difference_type offset,
    const WeightedClustersIterator<T, useSpans>& iter) noexcept {
  auto newIter = iter;
  newIter.index +=
      static_cast<
          typename WeightedClustersIterator<T, useSpans>::difference_type>(
          newIter.clusters) *
      offset;
  return newIter;
}

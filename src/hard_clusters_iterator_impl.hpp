#pragma once

#include "hard_cluster.hpp"
#include "hard_cluster_wrapper.hpp"
#include "hard_clusters_iterator.hpp"

#include <range/v3/algorithm.hpp>

template <typename Cluster, bool complemented>
HardClustersIterator<Cluster, complemented>::HardClustersIterator(
    Cluster &cluster, index_type offset) noexcept
    : wrapper(cluster, offset) {}

template <typename Cluster, bool complemented>
bool HardClustersIterator<Cluster, complemented>::operator<(
    const HardClustersIterator &other) const noexcept {
  assert(wrapper._cluster == other.wrapper._cluster);
  return wrapper._index < other.wrapper._index;
}

template <typename Cluster, bool complemented>
bool HardClustersIterator<Cluster, complemented>::operator==(
    const HardClustersIterator &other) const noexcept {
  assert(wrapper._cluster == other.wrapper._cluster);
  return wrapper._index == other.wrapper._index;
}

template <typename Cluster, bool complemented>
bool HardClustersIterator<Cluster, complemented>::operator!=(
    const HardClustersIterator &other) const noexcept {
  assert(wrapper._cluster == other.wrapper._cluster);
  return wrapper._index != other.wrapper._index;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator-(
    const HardClustersIterator &other) const noexcept -> difference_type {
  assert(wrapper._cluster == other.wrapper._cluster);
  return static_cast<difference_type>(wrapper._index) - other.wrapper._index;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator++() noexcept
    -> HardClustersIterator & {
  ++wrapper._index;
  return *this;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator++(int) noexcept
    -> HardClustersIterator {
  auto copy = *this;
  ++wrapper._index;
  return copy;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator--() noexcept
    -> HardClustersIterator & {
  --wrapper._index;
  return *this;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator--(int) noexcept
    -> HardClustersIterator {
  auto copy = *this;
  --wrapper._index;
  return copy;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator+=(
    difference_type offset) noexcept -> HardClustersIterator & {
  wrapper._index += offset;
  return *this;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator-=(
    difference_type offset) noexcept -> HardClustersIterator & {
  wrapper._index -= offset;
  return *this;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator+(
    difference_type offset) const noexcept -> HardClustersIterator {
  auto newIter = *this;
  newIter.wrapper._index += offset;
  return newIter;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator-(
    difference_type offset) const noexcept -> HardClustersIterator {
  auto newIter = *this;
  newIter.wrapper._index -= offset;
  return newIter;
}

template <typename _Cluster, bool _complemented>
HardClustersIterator<_Cluster, _complemented>
operator+(typename HardClustersIterator<_Cluster,
                                        _complemented>::difference_type offset,
          const HardClustersIterator<_Cluster, _complemented> &iter) noexcept {
  auto newIter = iter;
  newIter.wrapper._index += offset;
  return newIter;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator*() const noexcept
    -> reference {
  return wrapper;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator[](
    difference_type offset) const noexcept -> reference {
  auto newWrapper = wrapper;
  newWrapper._index += offset;
  return newWrapper;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator=(
    value_type const &rhs) const noexcept -> HardClustersIterator const & {
  assert(wrapper);
  assert(wrapper.size() == rhs.size());

  /*
  auto rhs_iter = std::begin(rhs);
  auto const end_rhs = std::end(rhs);
  auto wrapper_iter = wrapper.begin();

  for (; rhs_iter < end_rhs; ++rhs_iter, ++wrapper_iter) {
    auto element = *wrapper_iter;
    if (*rhs_iter)
      element.set();
    else
      element.clear();
  }
  asserT(wrapper_iter == wrapper.end());
  */

  ranges::copy(rhs, ranges::begin(*this));
  return this;
}

template <typename Cluster, bool complemented>
auto HardClustersIterator<Cluster, complemented>::operator=(
    value_type &&rhs) const noexcept -> HardClustersIterator const & {
  assert(wrapper);
  assert(wrapper.size() == rhs.size());

  ranges::move(rhs, ranges::begin(*this));
  return this;
}

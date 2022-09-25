#pragma once

#include "hard_cluster_element_iterator.hpp"
#include "hard_cluster_element_wrapper.hpp"
#include "hard_cluster_wrapper.hpp"

template <typename Cluster, bool complemented>
HardClusterElementIterator<Cluster, complemented>::HardClusterElementIterator(
    hard_cluster_wrapper_type const &cluster, std::size_t elementIndex) noexcept
    : wrapper(cluster, elementIndex) {}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator++() noexcept
    -> self & {
  ++wrapper._element;
  return *this;
}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator++(int) noexcept
    -> self {
  auto copy = *this;
  ++wrapper._element;
  return copy;
}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator--() noexcept
    -> self & {
  --wrapper._element;
  return *this;
}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator--(int) noexcept
    -> self {
  auto copy = *this;
  --wrapper._element;
  return copy;
}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator+=(
    difference_type offset) noexcept -> self & {
  wrapper._element += offset;
  return *this;
}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator-=(
    difference_type offset) noexcept -> self & {
  wrapper._element -= offset;
  return *this;
}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator+(
    difference_type offset) const noexcept -> self {
  auto newIter = *this;
  newIter._element += offset;
  return newIter;
}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator-(
    difference_type offset) const noexcept -> self {
  auto newIter = *this;
  newIter._element -= offset;
  return newIter;
}

template <typename _Cluster, bool _complemented>
HardClusterElementIterator<_Cluster, _complemented> operator+(
    typename HardClusterElementIterator<_Cluster,
                                        _complemented>::difference_type offset,
    const HardClusterElementIterator<_Cluster, _complemented> &iter) noexcept {
  auto newIter = iter;
  newIter._element += offset;
  return newIter;
}

template <typename Cluster, bool complemented>
bool HardClusterElementIterator<Cluster, complemented>::operator<(
    const self &other) const noexcept {
  return wrapper._element < other.wrapper._element;
}

template <typename Cluster, bool complemented>
bool HardClusterElementIterator<Cluster, complemented>::operator>(
    const self &other) const noexcept {
  return wrapper._element > other.wrapper._element;
}

template <typename Cluster, bool complemented>
bool HardClusterElementIterator<Cluster, complemented>::operator>=(
    const self &other) const noexcept {
  return wrapper._element >= other.wrapper._element;
}

template <typename Cluster, bool complemented>
bool HardClusterElementIterator<Cluster, complemented>::operator<=(
    const self &other) const noexcept {
  return wrapper._element <= other.wrapper._element;
}

template <typename Cluster, bool complemented>
bool HardClusterElementIterator<Cluster, complemented>::operator==(
    const self &other) const noexcept {
  return wrapper._element == other.wrapper._element;
}

template <typename Cluster, bool complemented>
bool HardClusterElementIterator<Cluster, complemented>::operator!=(
    const self &other) const noexcept {
  return wrapper._element != other.wrapper._element;
}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator-(
    const self &other) const noexcept -> difference_type {
  return wrapper._element - other.wrapper._element;
}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator*()
    const noexcept -> reference {
  return wrapper;
}

template <typename Cluster, bool complemented>
auto HardClusterElementIterator<Cluster, complemented>::operator[](
    difference_type offset) const noexcept -> reference {
  auto newWrapper = wrapper;
  newWrapper._element += offset;
  return newWrapper;
}

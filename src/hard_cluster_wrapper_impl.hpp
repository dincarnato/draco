#pragma once

#include "hard_cluster.hpp"
#include "hard_cluster_element_iterator.hpp"
#include "hard_cluster_element_wrapper.hpp"
#include "hard_cluster_wrapper.hpp"
#include "hard_clusters_iterator.hpp"

#include <cassert>
#include <limits>
#include <range/v3/algorithm.hpp>

template <typename Cluster, bool complemented>
HardClusterWrapper<Cluster, complemented>::HardClusterWrapper(
    Cluster &cluster, index_type index) noexcept
    : _cluster(&cluster), _index(index) {}

template <typename Cluster, bool complemented>
auto HardClusterWrapper<Cluster, complemented>::operator[](
    std::size_t index) const noexcept -> element_wrapper {
  assert(_index < _cluster->nClusters);
  assert(index < _cluster->_clusters.size());
  return {*this, index};
}

template <typename Cluster, bool complemented>
void HardClusterWrapper<Cluster, complemented>::set(
    std::size_t index) const noexcept {
  static_assert(not complemented);

  assert(_index < _cluster->nClusters);
  assert(index < _cluster->_clusters.size());
  _cluster->clusters[index] = _index;
}

template <typename Cluster, bool complemented>
void HardClusterWrapper<Cluster, complemented>::clear(
    std::size_t index) const noexcept {
  static_assert(not complemented);

  assert(_index < _cluster->nClusters);
  assert(index < _cluster->_clusters.size());
  _cluster->clusters[index] = std::numeric_limits<index_type>::max();
}

template <typename Cluster, bool complemented>
void HardClusterWrapper<Cluster, complemented>::moveTo(
    std::size_t index, std::size_t otherClusterIndex) const noexcept {
  static_assert(not complemented);

  assert(_index < _cluster->nClusters);
  assert(index < _cluster->_clusters.size());
  _cluster->clusters[index] = otherClusterIndex;
}

template <typename Cluster, bool complemented>
auto HardClusterWrapper<Cluster, complemented>::begin() const noexcept
    -> iterator {
  return {*this, 0};
}

template <typename Cluster, bool complemented>
auto HardClusterWrapper<Cluster, complemented>::end() const noexcept
    -> iterator {
  return {*this, _cluster->_clusters.size()};
}

template <typename Cluster, bool complemented>
auto HardClusterWrapper<Cluster, complemented>::index() const noexcept
    -> index_type {
  return _index;
}

template <typename Cluster, bool complemented>
std::size_t HardClusterWrapper<Cluster, complemented>::size() const noexcept {
  assert(_cluster);
  return _cluster->getElementsSize();
}

template <typename Cluster, bool complemented>
HardClusterWrapper<Cluster, complemented>::operator concrete_type() const
    noexcept(false) {
  assert(_cluster->nClusters > 0);
  /*
  auto const elements_size = _cluster->_clusters.size();
  std::vector<bool> out(elements_size, false);

  auto iter = std::begin(out);
  auto element_iter = begin();
  auto const end_iter = end();

  for (; element_iter < end_iter; ++element_iter, ++iter) {
    if (*element_iter)
      *iter = true;
  }

  return out;
  */
  return concrete_type(*this);
}

template <typename Cluster, bool complemented>
auto HardClusterWrapper<Cluster, complemented>::operator=(
    concrete_type const &rhs) const noexcept -> HardClusterWrapper const & {
  assert(_cluster->nClusters > 0);
  ranges::copy(rhs, ranges::begin(*this));

  return *this;
}

template <typename Cluster, bool complemented>
auto HardClusterWrapper<Cluster, complemented>::operator=(
    concrete_type &&rhs) const noexcept -> HardClusterWrapper const & {
  assert(_cluster->nClusters > 0);
  ranges::move(rhs, ranges::begin(*this));

  return *this;
}

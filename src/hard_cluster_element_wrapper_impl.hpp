#pragma once

#include "hard_cluster_element_wrapper.hpp"
#include "hard_cluster_wrapper.hpp"

#include <cassert>
#include <limits>

template <typename Cluster, bool complemented>
template <typename OtherCluster>
HardClusterElementWrapper<Cluster, complemented>::HardClusterElementWrapper(
    const HardClusterWrapper<OtherCluster, complemented>& clusterWrapper,
    std::size_t elementIndex) noexcept
    : _element(&clusterWrapper._cluster->_clusters[elementIndex]),
      clusterIndex(clusterWrapper._index)
#ifndef NDEBUG
      ,
      _cluster(clusterWrapper._cluster)
#endif
{
  static_assert(std::is_same_v<std::decay_t<OtherCluster>, cluster_type>);
}

template <typename Cluster, bool complemented>
auto
HardClusterElementWrapper<Cluster, complemented>::operator=(bool value) const
    noexcept -> HardClusterElementWrapper const& {
  if (value)
    set();
  else
    clear();
}

template <typename Cluster, bool complemented>
bool
HardClusterElementWrapper<Cluster, complemented>::get() const noexcept {
  assert(_element < _cluster->_clusters.data() + _cluster->_clusters.size());

  if constexpr (not complemented)
    return *_element == clusterIndex;
  else
    return *_element != clusterIndex;
}

template <typename Cluster, bool complemented>
void
HardClusterElementWrapper<Cluster, complemented>::set() const noexcept {
  static_assert(not complemented);

  assert(_element < _cluster->_clusters.data() + _cluster->_clusters.size());
  *_element = clusterIndex;
}

template <typename Cluster, bool complemented>
void
HardClusterElementWrapper<Cluster, complemented>::clear() const noexcept {
  static_assert(not complemented);

  assert(_element < _cluster->_clusters.data() + _cluster->_clusters.size());
  *_element = std::numeric_limits<index_type>::max();
}

template <typename Cluster, bool complemented>
HardClusterElementWrapper<Cluster, complemented>::operator bool() const
    noexcept {
  return get();
}

template <typename Cluster, bool complemented>
void
HardClusterElementWrapper<Cluster, complemented>::swap(
    HardClusterElementWrapper other) const noexcept {
  static_assert(not complemented);
  static_assert(not other.isComplemented);

  assert(_cluster == other._cluster);
  assert(_element < _cluster->_clusters.data() + _cluster->_clusters.size());
  assert(other._element <
         other._cluster->_clusters.data() + other._cluster->_clusters.size());

  if (_element != other._element) {
    if (get()) {
      if (other.get())
        std::swap(*_element, *other._element);
      else {
        if (clusterIndex != other.clusterIndex)
          *_element = other.clusterIndex;
        else
          clear();

        other.set();
      }
    } else {
      if (other.get()) {
        if (other.clusterIndex != clusterIndex)
          *other._element = clusterIndex;
        else
          other.clear();
      } else
        other.set();

      set();
    }
  } else {
    if (get())
      *_element = other.clusterIndex;
    else if (other.get())
      *_element = clusterIndex;
  }
}

template <typename Cluster, bool complemented>
void
HardClusterElementWrapper<Cluster, complemented>::swapMatched(
    HardClusterElementWrapper other) const noexcept {
  static_assert(not complemented);
  static_assert(not other.isComplemented);

  assert(_cluster == other._cluster);
  assert(_element < _cluster->_clusters.data() + _cluster->_clusters.size());
  assert(other._element <
         other._cluster->_clusters.data() + other._cluster->_clusters.size());
  assert(get());
  assert(other.get());

  if (_element != other._element)
    std::swap(*_element, *other._element);
}

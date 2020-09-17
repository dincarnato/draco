#pragma once

#include "hard_cluster_wrapper.hpp"
#include "hard_clusters_base_complement.hpp"
#include "hard_clusters_element_wrapper.hpp"
#include "hard_clusters_wrapper.hpp"

#include <cassert>

template <typename Base>
HardClustersBaseComplement<Base>::HardClustersBaseComplement(
    hard_clusters_base_type& clusters) noexcept
    : _clusters(&clusters) {}

template <typename Base>
auto
HardClustersBaseComplement<Base>::cluster(index_type index) const noexcept
    -> cluster_wrapper_type {
  assert(index < _clusters->nClusters);
  return {*_clusters, index};
}

template <typename Base>
auto HardClustersBaseComplement<Base>::operator[](std::size_t index) const
    noexcept -> element_wrapper {
  return {*_clusters, index};
}

template <typename Base>
auto
HardClustersBaseComplement<Base>::clusters() const noexcept
    -> clusters_wrapper_type {
  return {*_clusters};
}

template <typename Base>
auto
HardClustersBaseComplement<Base>::complement() const noexcept
    -> hard_clusters_base_type& {
  return *_clusters;
}

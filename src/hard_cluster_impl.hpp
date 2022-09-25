#pragma once

#include "hard_cluster.hpp"
#include "hard_cluster_wrapper.hpp"

#include <range/v3/algorithm.hpp>

template <typename Cluster, bool complemented>
inline HardCluster::HardCluster(
    HardClusterWrapper<Cluster, complemented> const &wrapper) noexcept(false) {
  init_from_cluser(wrapper);
}

template <typename Cluster, bool complemented>
inline HardCluster::HardCluster(
    HardClusterWrapper<Cluster, complemented> &&wrapper) noexcept(false) {
  init_from_cluser(std::move(wrapper));
}

template <typename Wrapper>
void HardCluster::init_from_wrapper(Wrapper &&wrapper) noexcept(false) {
  base_type::resize(wrapper.size());
  ranges::copy(wrapper, ranges::begin(*this));
}

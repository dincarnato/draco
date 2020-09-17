#pragma once

#include <vector>

template <typename, bool>
struct HardClusterWrapper;

struct HardCluster : std::vector<bool> {
  using base_type = std::vector<bool>;

  using base_type::base_type;

  template <typename Cluster, bool complemented>
  HardCluster(
      HardClusterWrapper<Cluster, complemented> const& wrapper) noexcept(false);

  template <typename Cluster, bool complemented>
  HardCluster(HardClusterWrapper<Cluster, complemented>&& wrapper) noexcept(
      false);

private:
  template <typename Wrapper>
  void init_from_wrapper(Wrapper&& wrapper) noexcept(false);
};

#include "hard_cluster_impl.hpp"

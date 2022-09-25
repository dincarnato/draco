#pragma once

#include "clusters_traits.hpp"
#include "weighted_clusters_cluster_wrapper.hpp"
#include "weighted_clusters_clusters_wrapper.hpp"
#include "weighted_clusters_iterator.hpp"
#include "weighted_clusters_span.hpp"

#include <vector>

class WeightedClusters {
  template <typename, bool> friend class WeightedClustersIterator;
  template <typename> friend class WeightedClustersClusterWrapper;
  template <typename> friend class WeightedClustersSpan;
  template <typename> friend class WeightedClustersClusterIterator;

  using weight_type = weighted_clusters_weight_type;
  using weights_type = weighted_clusters_weights_type;

public:
  using iterator = WeightedClustersIterator<weight_type>;
  using const_iterator = WeightedClustersIterator<const weight_type>;
  using span_type = WeightedClustersSpan<weight_type>;
  using const_span_type = WeightedClustersSpan<const weight_type>;
  using cluster_wrapper = WeightedClustersClusterWrapper<weight_type>;
  using const_cluster_wrapper =
      WeightedClustersClusterWrapper<const weight_type>;
  using clusters_wrapper = WeightedClustersClustersWrapper<weight_type>;
  using const_clusters_wrapper =
      WeightedClustersClustersWrapper<const weight_type>;

  WeightedClusters() = default;
  WeightedClusters(std::size_t nElements, std::size_t nClusters,
                   bool setDefaultCluster = true);

  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  span_type operator[](std::size_t index);
  const_span_type operator[](std::size_t index) const;

  cluster_wrapper cluster(std::size_t index);
  const_cluster_wrapper cluster(std::size_t index) const;

  clusters_wrapper clusters();
  const_clusters_wrapper clusters() const;

  void setOnlyFirstCluster();
  void removeWeights();

  std::size_t getClustersSize() const;
  std::size_t getElementsSize() const;

  WeightedClusters complement() const;

  bool operator==(const WeightedClusters &other) const noexcept;
  bool operator<(const WeightedClusters &other) const noexcept;

private:
  std::size_t elements = 0;
  std::size_t _clusters = 0;
  weights_type weights;
};

static_assert(is_clusters_v<WeightedClusters>);

#include "weighted_clusters_impl.hpp"

#pragma once

#include "weighted_clusters.hpp"

#include <algorithm>
#include <ranges>
#include <stdexcept>

inline WeightedClusters::WeightedClusters(std::size_t nElements,
                                          std::size_t nClusters,
                                          bool setDefaultCluster)
    : elements(nElements), _clusters(nClusters),
      weights(nElements * nClusters, 0.f) {
  if (setDefaultCluster) {
    for (auto iter = std::begin(weights); iter < std::end(weights);
         iter += static_cast<typename weights_type::difference_type>(_clusters))
      *iter = 1.f;
  }
}

inline WeightedClusters::WeightedClusters(
    std::initializer_list<std::initializer_list<WeightedClusters::weight_type>>
        weights)
    : elements(std::empty(weights) ? 0
                                   : std::size(*std::ranges::begin(weights))),
      _clusters(std::size(weights)), weights(elements * _clusters) {
  if (std::ranges::any_of(
          weights | std::views::drop(1),
          [&](auto cluster_size) { return cluster_size != elements; },
          [](auto &&cluster_weights) { return std::size(cluster_weights); })) {
    throw std::runtime_error("invalid initializer list for WeightedClusters");
  }

  std::ranges::for_each(
      std::views::zip(std::views::iota(0uz) |
                          std::views::transform([this](auto cluster_index) {
                            return this->weights |
                                   std::views::drop(cluster_index) |
                                   std::views::stride(_clusters);
                          }),
                      weights),
      [](auto &&tuple) {
        auto &&[cluster_to_assign, cluster] = tuple;
        std::ranges::copy(cluster, std::ranges::begin(cluster_to_assign));
      });
}

inline auto WeightedClusters::begin() -> iterator { return iterator{*this, 0}; }

inline auto WeightedClusters::end() -> iterator {
  return iterator(*this, static_cast<std::ptrdiff_t>(elements));
}

inline auto WeightedClusters::begin() const -> const_iterator {
  return const_iterator{*this, 0};
}

inline auto WeightedClusters::end() const -> const_iterator {
  return const_iterator(*this, static_cast<std::ptrdiff_t>(elements));
}

inline auto WeightedClusters::operator[](std::size_t index) -> span_type {
  assert(index < elements);
  return {weights, _clusters,
          static_cast<std::ptrdiff_t>(index) *
              static_cast<std::ptrdiff_t>(_clusters)};
}

inline auto WeightedClusters::operator[](std::size_t index) const
    -> const_span_type {
  assert(index < elements);
  return {weights, _clusters,
          static_cast<std::ptrdiff_t>(index) *
              static_cast<std::ptrdiff_t>(_clusters)};
}

inline auto WeightedClusters::cluster(std::size_t index) -> cluster_wrapper {
  return {*this, static_cast<std::ptrdiff_t>(index)};
}

inline auto WeightedClusters::cluster(std::size_t index) const
    -> const_cluster_wrapper {
  return {*this, static_cast<std::ptrdiff_t>(index)};
}

inline void WeightedClusters::setOnlyFirstCluster() {
  std::fill(std::begin(weights), std::end(weights), 0);
  auto endWeights = std::end(weights);
  for (auto iter = std::begin(weights); iter < endWeights;
       iter += static_cast<typename weights_type::difference_type>(_clusters))
    *iter = 1.f;
}

inline void WeightedClusters::removeWeights() {
  std::fill(std::begin(weights), std::end(weights), 0);
}

inline auto WeightedClusters::clusters() -> clusters_wrapper {
  return clusters_wrapper{*this};
}

inline auto WeightedClusters::clusters() const -> const_clusters_wrapper {
  return const_clusters_wrapper{*this};
}

inline std::size_t WeightedClusters::getClustersSize() const {
  return _clusters;
}

inline std::size_t WeightedClusters::getElementsSize() const {
  return elements;
}

inline WeightedClusters WeightedClusters::complement() const {
  WeightedClusters complement(elements, _clusters, false);
  std::transform(std::begin(weights), std::end(weights),
                 std::begin(complement.weights),
                 [](auto weight) { return 1.f - weight; });

  return complement;
}

inline bool
WeightedClusters::operator==(const WeightedClusters &other) const noexcept {
  return elements == other.elements and _clusters == other._clusters and
         std::equal(std::begin(weights), std::end(weights),
                    std::begin(other.weights), std::end(other.weights));
}

inline bool
WeightedClusters::operator<(const WeightedClusters &other) const noexcept {
  if (elements < other.elements)
    return true;
  else if (elements == other.elements) {
    if (_clusters < other._clusters)
      return true;
    else if (_clusters == other._clusters)
      return weights < other.weights;
  }

  return false;
}

static_assert(std::random_access_iterator<typename WeightedClusters::iterator>);
static_assert(
    std::random_access_iterator<typename WeightedClusters::const_iterator>);

static_assert(std::random_access_iterator<
              typename WeightedClusters::cluster_wrapper::iterator>);
static_assert(std::random_access_iterator<
              typename WeightedClusters::const_cluster_wrapper::iterator>);

static_assert(std::random_access_iterator<
              typename WeightedClusters::clusters_wrapper::iterator>);
static_assert(std::random_access_iterator<
              typename WeightedClusters::const_clusters_wrapper::iterator>);

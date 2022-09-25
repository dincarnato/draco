#pragma once

#include <type_traits>
#include <vector>

template <typename T> class WeightedClustersSpan;

template <typename T> struct WeightedClustersConcreteSpan : std::vector<T> {
  using base_type = std::vector<T>;

  WeightedClustersConcreteSpan() = default;
  template <typename U>
  WeightedClustersConcreteSpan(WeightedClustersSpan<U> const &span) noexcept(
      false);
  template <typename U>
  WeightedClustersConcreteSpan(WeightedClustersSpan<U> &&span) noexcept(false);

private:
  template <typename Span> void init_from_span(Span &&span) noexcept(false);
};

template <typename T>
WeightedClustersConcreteSpan(WeightedClustersSpan<T> const &)
    -> WeightedClustersConcreteSpan<std::decay_t<T>>;

template <typename T>
WeightedClustersConcreteSpan(WeightedClustersSpan<T> &&)
    -> WeightedClustersConcreteSpan<std::decay_t<T>>;

#include "weighted_clusters_concrete_span_impl.hpp"

#pragma once

#include "weighted_clusters_concrete_span.hpp"
#include "weighted_clusters_span.hpp"

#include <range/v3/algorithm.hpp>

template <typename T>
template <typename U>
WeightedClustersConcreteSpan<T>::WeightedClustersConcreteSpan(
    WeightedClustersSpan<U> const &span) noexcept(false) {
  static_assert(std::is_same_v<std::decay_t<U>, T>);
  init_from_span(span);
}

template <typename T>
template <typename U>
WeightedClustersConcreteSpan<T>::WeightedClustersConcreteSpan(
    WeightedClustersSpan<U> &&span) noexcept(false) {
  static_assert(std::is_same_v<std::decay_t<U>, T>);
  init_from_span(std::move(span));
}

template <typename T>
template <typename Span>
void WeightedClustersConcreteSpan<T>::init_from_span(Span &&span) noexcept(
    false) {
  base_type::resize(span.size());
  ranges::copy(span, std::begin(*this));
}

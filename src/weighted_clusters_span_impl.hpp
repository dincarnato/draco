#pragma once

#include "weighted_clusters_span.hpp"

#include <algorithm>
#include <cassert>
#include <limits>

template <typename T>
WeightedClustersSpan<T>::WeightedClustersSpan(weights_reference weights,
                                              std::size_t span,
                                              std::ptrdiff_t index) noexcept
    : _data(weights.data()), index(index), _size(weights.size()), _span(span) {
  assert(_data);
}

template <typename T>
auto WeightedClustersSpan<T>::operator=(concrete_type const &concrete_span)
    const noexcept -> WeightedClustersSpan const & {
  assert(concrete_span.size() == _span);
  std::ranges::copy(concrete_span, std::ranges::begin(*this));
  return *this;
}

template <typename T>
auto WeightedClustersSpan<T>::operator=(concrete_type &&concrete_span)
    const noexcept -> WeightedClustersSpan const & {
  assert(concrete_span.size() == _span);
  std::ranges::move(concrete_span, std::ranges::begin(*this));
  return *this;
}

template <typename T>
WeightedClustersSpan<T>::operator concrete_type() const noexcept(false) {
  return concrete_type(*this);
}

template <typename T>
auto WeightedClustersSpan<T>::operator[](std::size_t index) const noexcept
    -> reference {
  assert(_data);
  return _data[getIndex(index)];
}

template <typename T>
auto WeightedClustersSpan<T>::begin() const noexcept -> iterator {
  assert(_data);
  return _data + getIndex();
}

template <typename T>
auto WeightedClustersSpan<T>::end() const noexcept -> iterator {
  assert(_data);
  return _data + getIndex(_span);
}

template <typename T>
std::size_t WeightedClustersSpan<T>::data_size() const noexcept {
  assert(_data);
  return _size + getIndex();
}

template <typename T>
std::ptrdiff_t
WeightedClustersSpan<T>::getIndex(std::ptrdiff_t offset) const noexcept {
  return std::min(std::max(static_cast<std::ptrdiff_t>(0), index + offset),
                  static_cast<std::ptrdiff_t>(_size));
}

template <typename T>
std::ptrdiff_t
WeightedClustersSpan<T>::getIndex(std::size_t offset) const noexcept {
  assert(offset <= std::numeric_limits<std::ptrdiff_t>::max());
  return getIndex(static_cast<std::ptrdiff_t>(offset));
}

template <typename T> const T *WeightedClustersSpan<T>::data() const noexcept {
  return data;
}

template <typename T>
auto WeightedClustersSpan<T>::operator*() const noexcept -> reference {
  assert(_data);
  return _data[getIndex()];
}

template <typename T>
std::size_t WeightedClustersSpan<T>::span_size() const noexcept {
  return _span;
}

static_assert(
    std::random_access_iterator<WeightedClustersSpan<float>::iterator>);
static_assert(
    std::random_access_iterator<WeightedClustersSpan<const float>::iterator>);

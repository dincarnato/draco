#pragma once

#include "clusters_traits.hpp"

#include <cstddef>
#include <type_traits>

class WeightedClusters;
template <typename T>
struct WeightedClustersConcreteSpan;

template <typename T>
class WeightedClustersSpan {
  using weights_reference =
      std::conditional_t<std::is_const_v<T>,
                         const weighted_clusters_weights_type&,
                         weighted_clusters_weights_type&>;
  static_assert(
      std::is_same_v<weighted_clusters_weight_type, std::decay_t<T>>,
      "weighted_clusters_weight_type type must be the same as decayed T");
  using concrete_type = WeightedClustersConcreteSpan<std::decay_t<T>>;

public:
  using iterator = T*;
  using reference = std::conditional_t<std::is_rvalue_reference_v<T>, T, T&>;

  explicit WeightedClustersSpan() = default;
  WeightedClustersSpan(weights_reference weights, std::size_t span,
                       std::ptrdiff_t index) noexcept;

  WeightedClustersSpan const&
  operator=(concrete_type const& concrete_span) const noexcept;
  WeightedClustersSpan const& operator=(concrete_type&& concrete_span) const
      noexcept;

  operator concrete_type() const noexcept(false);

  reference operator[](std::size_t index) const noexcept;
  reference operator*() const noexcept;

  iterator begin() const noexcept;
  iterator end() const noexcept;

  std::size_t data_size() const noexcept;
  std::size_t span_size() const noexcept;
  const T* data() const noexcept;

private:
  T* _data;

public:
  std::ptrdiff_t index;

private:
  std::size_t _size;
  std::size_t _span;

  inline std::ptrdiff_t getIndex(std::ptrdiff_t offset = 0) const noexcept;
  inline std::ptrdiff_t getIndex(std::size_t offset) const noexcept;
};

#include "weighted_clusters_span_impl.hpp"

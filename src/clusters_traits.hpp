#pragma once

#include <iterator>
#include <range/v3/iterator/concepts.hpp>
#include <type_traits>
#include <vector>

using weighted_clusters_weight_type = float;
using weighted_clusters_weights_type =
    std::vector<weighted_clusters_weight_type>;

template <typename> struct HardClustersBase;

template <typename, typename = void> struct HardClustersTraits {};

template <typename Type> struct HardClustersTraits<HardClustersBase<Type>> {
  using cluster_type = HardClustersBase<Type>;
  using index_type = Type;
};

template <typename Type>
struct HardClustersTraits<const HardClustersBase<Type>> {
  using cluster_type = HardClustersBase<Type>;
  using index_type = Type;
};

template <typename T, typename = void>
struct is_cluster_wrapper : std::bool_constant<false> {};

template <typename T>
struct is_cluster_wrapper<
    T, std::enable_if_t<ranges::RandomAccessIterator<typename T::iterator>>>
    : std::bool_constant<true> {};

template <typename T>
constexpr bool is_cluster_wrapper_v = is_cluster_wrapper<T>::value;

template <typename T, typename = void>
struct is_clusters : std::bool_constant<false> {};

template <typename T>
struct is_clusters<T, std::void_t<decltype(std::declval<T>().clusters()),
                                  decltype(std::declval<T>().cluster(0)),
                                  decltype(std::declval<T>().complement()),
                                  std::enable_if_t<is_cluster_wrapper_v<
                                      decltype(std::declval<T>().cluster(0))>>>>
    : std::bool_constant<true> {};

template <typename T> constexpr bool is_clusters_v = is_clusters<T>::value;

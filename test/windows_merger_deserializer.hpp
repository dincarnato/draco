#pragma once

#include "json_deserializer.hpp"
#include "weighted_clusters.hpp"
#include "windows_merger_window.hpp"
#include <boost/callable_traits.hpp>

#include <cmath>
#include <range/v3/core.hpp>

namespace ct = boost::callable_traits;

#define REFLECT_IF_VISIT(name, value, variable, visitor)                       \
  if (name == #variable)                                                       \
    std::visit(                                                                \
        [&](auto &&v) { visitor(variable, std::forward<decltype(v)>(v)); },    \
        value);

#define REFLECT_ELSE_IF_VISIT(name, value, variable, visitor)                  \
  else REFLECT_IF_VISIT(name, value, variable, visitor)

#define REFLECT_IF_FUN(name, value, variable, fun)                             \
  if (name == #variable)                                                       \
    fun(variable,                                                              \
        std::get<                                                              \
            std::decay_t<std::tuple_element_t<1, ct::args_t<decltype(fun)>>>>( \
            value));

#define REFLECT_ELSE_IF_FUN(name, value, variable, fun)                        \
  else REFLECT_IF_FUN(name, value, variable, fun)

#define REFLECT_IF(name, value, variable)                                      \
  if (name == #variable)                                                       \
  std::visit([&](auto &&value) { variable = value }, value)

#define REFLECT_ELSE_IF(name, value, variable)                                 \
  else REFLECT_IF(name, value, variable)

#define REFLECT_IF_MOVE(name, value, variable)                                 \
  if (name == #variable)                                                       \
  std::visit([&](auto &&value) { variable = std::move(value) }, value)

#define REFLECT_ELSE_IF_MOVE(name, value, variable)                            \
  else REFLECT_IF_MOVE(name, value, variable)

void weighted_clusters_from_json(WeightedClusters &weighted_clusters,
                                 json_deserializer::Array const &json_data) {
  namespace jsonde = json_deserializer;
  assert(std::size(json_data) == weighted_clusters.getClustersSize());

  auto weighted_clusters_clusters = weighted_clusters.clusters();
  auto weighted_cluster_iter = weighted_clusters_clusters.begin();
  auto const weighted_clusters_end = weighted_clusters_clusters.end();
  auto json_cluster_iter = ranges::begin(json_data);

  for (; weighted_cluster_iter < weighted_clusters_end;
       ++weighted_cluster_iter, ++json_cluster_iter) {
    auto &&weighted_cluster = *weighted_cluster_iter;
    auto &&json_cluster = std::get<jsonde::Array>(*json_cluster_iter);

    ranges::transform(
        json_cluster, weighted_cluster.begin(), [](auto &&json_weight) {
          return static_cast<float>(
              std::get<double>(std::get<jsonde::Number>(json_weight)));
        });
  }
}

void coverages_from_json(std::vector<unsigned> &coverages,
                         json_deserializer::Array const &json_data) noexcept {
  coverages.resize(json_data.size());
  ranges::transform(
      json_data, ranges::begin(coverages), [](auto &&json_number) {
        return static_cast<unsigned>(
            std::get<long>(std::get<json_deserializer::Number>(json_number)));
      });
}

template <typename T>
void from_json_number(T &t, json_deserializer::Number const &json_number) {
  std::visit([&](auto &&value) { t = static_cast<T>(value); }, json_number);
}

template <typename T>
void from_json_number_value(T &t, json_deserializer::Value const &json_value) {
  namespace jsonde = json_deserializer;

  if (auto json_number_ptr = std::get_if<jsonde::Number>(&json_value);
      json_number_ptr != nullptr)
    std::visit([&](auto &&value) { t = static_cast<T>(value); },
               *json_number_ptr);
  else if (std::holds_alternative<jsonde::Null>(json_value)) {
    if constexpr (std::is_floating_point_v<T>)
      t = static_cast<T>(std::nan(""));
    else
      t = std::numeric_limits<T>::max();
  } else
    throw std::invalid_argument("cannot parse a non-number type as a number");
}

namespace json_deserializer {

bool deserialize(std::istream &is,
                 std::tuple<WeightedClusters, std::vector<unsigned>,
                            unsigned short> &window) {
  Object raw_obj;
  deserialize(is, raw_obj);

  auto &[weighted_clusters, coverages, start_base_index] = window;

  for (auto &&[name, value] : raw_obj) {
    REFLECT_IF_FUN(name, value, weighted_clusters, weighted_clusters_from_json)
    REFLECT_ELSE_IF_FUN(name, value, coverages, coverages_from_json)
    REFLECT_ELSE_IF_FUN(name, value, start_base_index,
                        from_json_number<unsigned short>)
    else throw std::invalid_argument("cannot serialize data");
  }

  return true;
}

} // namespace json_deserializer

static std::tuple<WeightedClusters, std::vector<unsigned>, unsigned short>
deserialize_initial_window(std::istream &is, unsigned n_clusters,
                           unsigned n_elements) {
  namespace jsonde = json_deserializer;

  auto out = std::tuple(WeightedClusters(n_elements, n_clusters, false),
                        std::vector<unsigned>(n_elements),
                        static_cast<unsigned short>(0));
  jsonde::deserialize(is, out);
  return out;
}

namespace {

using Weights = std::vector<double>;
struct Window {
  std::vector<Weights> clusters_weights;
  std::vector<unsigned> coverages;
  unsigned short begin{};
  unsigned short end{};

  Window() = default;
  Window(unsigned char n_clusters)
      : clusters_weights(std::vector<Weights>(n_clusters)) {}

  void from_json(json_deserializer::Object const &obj) {
    auto &clusters = clusters_weights;

    for (auto &&[name, value] : obj) {
      REFLECT_IF_FUN(name, value, clusters, Window::clusters_from_json)
      REFLECT_ELSE_IF_FUN(name, value, coverages, coverages_from_json)
      REFLECT_ELSE_IF_FUN(name, value, begin, from_json_number<unsigned short>)
      REFLECT_ELSE_IF_FUN(name, value, end, from_json_number<unsigned short>)
      else throw std::invalid_argument("cannot serialize data");
    }
  }

private:
  static void clusters_from_json(std::vector<Weights> &clusters,
                                 json_deserializer::Array const &raw_clusters) {
    namespace jsonde = json_deserializer;
    assert(std::size(clusters) == std::size(raw_clusters));

    auto cluster_iter = ranges::begin(clusters);
    auto const clusters_end = ranges::end(clusters);
    auto raw_cluster_iter = ranges::begin(raw_clusters);
    for (; cluster_iter < clusters_end; ++cluster_iter, ++raw_cluster_iter) {
      auto &&cluster = *cluster_iter;
      auto &&raw_cluster = std::get<jsonde::Array>(*raw_cluster_iter);

      cluster.resize(raw_cluster.size());
      ranges::transform(
          raw_cluster, ranges::begin(cluster), [](auto &&json_number) {
            double weight;
            from_json_number(weight, std::get<jsonde::Number>(json_number));
            return weight;
          });
    }
  }
};

inline bool operator==(Window const &a,
                       windows_merger::WindowsMergerWindow const &b) noexcept {
  using clusters_size_type =
      windows_merger::WindowsMergerTraits::clusters_size_type;
  using bases_size_type = windows_merger::WindowsMergerTraits::bases_size_type;

  if (a.begin != b.begin_index() or a.end != b.end_index() or
      not ranges::equal(a.coverages, b.coverages()))
    return false;

  auto const clusters_size = a.clusters_weights.size();
  if (clusters_size != b.clusters_size())
    return false;

  if (a.end - a.begin == 0)
    return true;

  auto const bases_size = b.size();
  assert(bases_size == a.coverages.size());
  for (clusters_size_type cluster_index = 0; cluster_index < clusters_size;
       ++cluster_index) {
    auto &&a_cluster_weights = a.clusters_weights[cluster_index];
    for (bases_size_type base_index = 0; base_index < bases_size;
         ++base_index) {
      TinyFraction const a_weight(a_cluster_weights[base_index]);
      if (a_weight != b[base_index].weight(cluster_index))
        return false;
    }
  }

  return true;
}

inline bool operator==(windows_merger::WindowsMergerWindow const &a,
                       Window const &b) noexcept {
  return operator==(b, a);
}

inline bool operator!=(Window const &a,
                       windows_merger::WindowsMergerWindow const &b) noexcept {
  return not(a == b);
}

inline bool operator!=(windows_merger::WindowsMergerWindow const &a,
                       Window const &b) noexcept {
  return not(a == b);
}

} // namespace

static std::vector<Window> deserialize_windows(std::istream &is,
                                               unsigned char n_clusters) {
  namespace jsonde = json_deserializer;

  jsonde::Array raw_windows;
  jsonde::deserialize(is, raw_windows);

  std::vector<Window> windows;
  windows.resize(raw_windows.size(), Window(n_clusters));

  auto raw_window_iter = ranges::cbegin(raw_windows);
  auto const raw_windows_end = ranges::cend(raw_windows);
  auto window_iter = ranges::begin(windows);

  for (; raw_window_iter < raw_windows_end; ++raw_window_iter, ++window_iter)
    window_iter->from_json(std::get<jsonde::Object>(*raw_window_iter));

  return windows;
}

static Window deserialize_window(std::istream &is, unsigned char n_clusters) {
  namespace jsonde = json_deserializer;

  jsonde::Object raw_window;
  jsonde::deserialize(is, raw_window);

  Window window(n_clusters);
  window.from_json(raw_window);

  return window;
}

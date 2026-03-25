#include "weighted_clusters.hpp"
#include "ringmap_data.hpp"

#include <algorithm>
#include <optional>

std::optional<WeightedClusters>
WeightedClusters::create_reduced(RingmapData const &filtered_ringmap) const {
  if (not filtered_ringmap.bases_filtered()) {
    return std::nullopt;
  }

  WeightedClusters reduced_weights(filtered_ringmap.data().cols_size(),
                                   _clusters, false);
  for (const auto &old_and_new_col : filtered_ringmap.old_cols_to_new()) {
    assert(old_and_new_col.first < elements);
    assert(old_and_new_col.second < reduced_weights.getElementsSize());

    auto &&element_weights = (*this)[old_and_new_col.first];
    auto &&reduced_element_weights = reduced_weights[old_and_new_col.second];
    assert(element_weights.span_size() == reduced_element_weights.span_size());

    std::ranges::copy(element_weights,
                      std::ranges::begin(reduced_element_weights));
  }

  return reduced_weights;
}

void WeightedClusters::copy_from_reduced(
    WeightedClusters const &reduced_weights,
    RingmapData const &filtered_ringmap) noexcept {
  for (const auto &old_and_new_col : filtered_ringmap.old_cols_to_new()) {
    assert(old_and_new_col.first < elements);
    assert(old_and_new_col.second < reduced_weights.getElementsSize());

    auto &&element_weights = (*this)[old_and_new_col.first];
    auto &&reduced_element_weights = reduced_weights[old_and_new_col.second];
    assert(element_weights.span_size() == reduced_element_weights.span_size());

    std::ranges::copy(reduced_element_weights,
                      std::ranges::begin(element_weights));
  }
}

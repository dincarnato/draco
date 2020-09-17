#pragma once

#include "windows_merger.hpp"

namespace windows_merger {

template <typename Weights, typename Coverages>
void
WindowsMerger::add_window(bases_size_type start_offset, Weights&& weights,
                          Coverages&& coverages) {
  static_assert(std::is_same_v<std::decay_t<Weights>, input_weights_type>);
  static_assert(std::is_same_v<std::decay_t<Coverages>, input_coverages_type>);

  using namespace exceptions;
  if (weights.getElementsSize() != std::size(coverages))
    throw InvalidWindowSize("weights and coverages have different sizes");

  if (weights.getClustersSize() != windows.clusters_size())
    throw InvalidClustersSize("weights object have a wrong number of clusters");

  if (bool expected = false; dequeueing.compare_exchange_strong(
          expected, true, std::memory_order::memory_order_acq_rel)) {
    transform_window_to_native(start_offset, weights, coverages);
    process_queue();
  } else {
    queue.emplace(start_offset, std::forward<Weights>(weights),
                  std::forward<Coverages>(coverages));
  }
}

} // namespace windows_merger

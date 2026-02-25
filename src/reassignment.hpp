#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

struct RingmapData;
struct PtbaOnReplicate;
struct Args;

namespace results {
struct Transcript;
} // namespace results

struct Reassignment {
  std::span<const std::vector<RingmapData>> replicates_splitted_ringmaps;
  std::span<RingmapData> filtered_ringmaps;
  std::span<const PtbaOnReplicate> ptba_on_replicate_results;
  std::span<std::optional<unsigned>> windows_max_clusters_constraints;
  results::Transcript *transcript_result;
  bool *stop;
  Args const *args;
  std::size_t window_index;
  unsigned window_size;
  bool allow_empty_patterns;

  void reassign_reads_with_weights() const;
};

#pragma once

#include "weighted_clusters.hpp"

#include <mutex>
#include <vector>

struct Window {
  unsigned short start_base;
  WeightedClusters weights;
  std::vector<unsigned> coverages;
};

class MutationMapTranscript;
class RingmapData;

namespace results {
struct Analysis;
struct Transcript;
} // namespace results

struct Args;

void handle_transcript(MutationMapTranscript const &transcript,
                       RingmapData &ringmapData,
                       results::Analysis &analysisResult, Args const &args,
                       std::optional<std::ofstream> &raw_n_clusters_stream,
                       std::mutex &raw_n_clusters_stream_mutex);

void merge_windows_and_add_window_results(
    std::vector<Window> const &windows,
    std::vector<std::vector<unsigned>> const &windows_reads_indices,
    RingmapData &ringmap_data, results::Transcript &transcript_result,
    Args const &args);

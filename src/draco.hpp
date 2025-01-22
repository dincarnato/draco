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
} // namespace results

struct Args;

void handle_transcript(MutationMapTranscript const &transcript,
                       RingmapData &ringmapData,
                       results::Analysis &analysisResult, Args const &args,
                       std::optional<std::ofstream> &raw_n_clusters_stream,
                       std::mutex &raw_n_clusters_stream_mutex);

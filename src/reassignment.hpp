#pragma once

#include "ringmap_data.hpp"

#include <cstddef>
#include <optional>
#include <random>
#include <span>

struct PtbaOnReplicate;
struct Args;

namespace results {
struct Window;
struct Transcript;
} // namespace results

namespace reassignment {
struct ReweightAndReassignWithExpectationMaximizationIteration;
using FractionResult = std::tuple<RingmapData::clusters_fraction_type,
                                  RingmapData::clusters_pattern_type,
                                  RingmapData::clusters_assignment_type>;
} // namespace reassignment

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
  void reweight_and_reassign_with_expectation_maximization() const;

protected:
  reassignment::FractionResult
  reweight_and_reassign_with_expectation_maximization_iteration(
      reassignment::ReweightAndReassignWithExpectationMaximizationIteration
          args) const;
};

namespace reassignment {

struct HandleFractionedReads {
  results::Window *window;
  RingmapData const *filtered_ringmap;
  RingmapData const *ringmap;
  bool *stop;
  results::Transcript *transcript_result;
  Args const *args;
  std::size_t replicate_index;
  std::size_t window_index;
  std::size_t window_size;
  std::span<const PtbaOnReplicate> ptba_on_replicate_results;
  std::span<std::optional<unsigned>> windows_max_clusters_constraints;
  FractionResult fractions_result;
  bool allow_empty_patterns;

  void operator()();
};

struct ReweightAndReassignWithExpectationMaximizationIteration {
  RingmapData const *filtered_ringmap;
  results::Window *window;
  std::mt19937 *rng;
  std::vector<std::uint32_t> *assignments_per_cluster;
  std::vector<double> *buffer;
  std::vector<std::uint32_t> *mapped_rows;
  std::vector<std::uint32_t> *clusters_reads_count;
  std::size_t replicate_index;
};

} // namespace reassignment

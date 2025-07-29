#pragma once

#include <algorithm>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "../window_clusters_with_confidence.hpp"
#include "jsonify_base.hpp"
#include "window.hpp"

namespace results {

struct WindowRange {
  std::size_t window_index_begin;
  std::size_t window_index_end;
};

struct Transcript {
  explicit constexpr Transcript(std::size_t replicates)
      : reads(replicates), coverages(replicates), windows(replicates),
        errors(replicates) {}

  std::string name;
  std::string sequence;
  std::vector<unsigned> reads;
  std::vector<std::optional<std::vector<unsigned>>> coverages;
  std::vector<std::optional<std::vector<Window>>> windows;
  std::vector<WindowRange> window_ranges;
  std::vector<WindowClustersWithConfidence> detected_clusters_with_confidence;
  std::vector<std::optional<std::string>> errors;
};

} // namespace results

namespace results {

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &jsonify_window_clusters_with_confidence(
    std::basic_ostream<CharT, Traits> &os,
    std::vector<WindowClustersWithConfidence> const
        &detected_clusters_with_confidence) {
  os << '[';
  auto write_data = [&](auto const &data) {
    os << '{';
    jsonify(os, "nClusters", data.n_clusters, "confidence", data.confidence,
            "start", data.start_base, "end", data.end_base);
    os << '}';
  };
  write_data(detected_clusters_with_confidence[0]);
  std::ranges::for_each(detected_clusters_with_confidence | std::views::drop(1),
                        [&](auto const &data) {
                          os << ',';
                          write_data(data);
                        });
  os << ']';

  return os;
}

template <typename CharT, typename Traits, typename T>
std::enable_if_t<std::is_convertible_v<std::decay_t<T>, Transcript>,
                 std::basic_ostream<CharT, Traits> &>
jsonify(std::basic_ostream<CharT, Traits> &os, T &&transcript) {
  os << '{';
  jsonify(os, "id", transcript.name, "sequence", transcript.sequence, "nReads",
          transcript.reads);

  if (std::ranges::any_of(transcript.coverages, [](auto const &maybeCoverages) {
        return maybeCoverages.has_value();
      })) {
    os << ',';
    jsonify(os, "preCoverages", transcript.coverages);
  }

  if (std::ranges::any_of(transcript.windows, [](auto const &maybeWindows) {
        return maybeWindows.has_value();
      })) {
    os << ',';
    jsonify(os, "windows", transcript.windows);
  }

  if (not transcript.detected_clusters_with_confidence.empty()) {
    os << ",\"replicateClusters\":";
    jsonify_window_clusters_with_confidence(
        os, transcript.detected_clusters_with_confidence);
  }

  if (std::ranges::any_of(transcript.errors, [](auto const &maybeError) {
        return maybeError.has_value();
      })) {
    os << ',';
    jsonify(os, "errors", transcript.errors);
  }

  return os << '}';
}

} // namespace results

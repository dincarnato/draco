#pragma once

#include <algorithm>
#include <optional>
#include <span>
#include <string>
#include <vector>

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
  std::vector<std::optional<std::string>> errors;
};

} // namespace results

namespace results {

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

  if (std::ranges::any_of(transcript.errors, [](auto const &maybeError) {
        return maybeError.has_value();
      })) {
    os << ',';
    jsonify(os, "errors", transcript.errors);
  }

  return os << '}';
}

} // namespace results

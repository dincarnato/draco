#pragma once

#include <optional>
#include <string>
#include <vector>

#include "jsonify_base.hpp"

namespace results {

struct Window;

struct Transcript {
  std::string name;
  std::string sequence;
  unsigned reads;
  std::optional<std::vector<unsigned>> coverages;
  std::optional<std::vector<Window>> windows;
  std::optional<std::string> error;
};

} // namespace results

namespace results {

template <typename CharT, typename Traits, typename T>
std::enable_if_t<std::is_convertible_v<std::decay_t<T>, Transcript>,
                 std::basic_ostream<CharT, Traits>&>
jsonify(std::basic_ostream<CharT, Traits>& os, T&& transcript) {
  os << '{';
  jsonify(os, "id", transcript.name, "sequence", transcript.sequence, "nReads",
          transcript.reads);

  if (transcript.coverages) {
    os << ',';
    jsonify(os, "preCoverage", transcript.coverages);
  }

  if (transcript.windows) {
    os << ',';
    jsonify(os, "windows", transcript.windows);
  }

  if (transcript.error) {
    os << ',';
    jsonify(os, "error", transcript.error);
  }

  return os << '}';
}

} // namespace results

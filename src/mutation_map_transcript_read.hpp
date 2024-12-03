#pragma once

#include <limits>
#include <vector>

struct MutationMapTranscriptRead {
  unsigned begin;
  unsigned end;
  std::vector<unsigned> indices;

  constexpr unsigned original_begin_index() const { return begin; }
  constexpr unsigned original_end_index() const { return end; }
  constexpr unsigned window_begin_index() const { return 0; }
  constexpr unsigned window_end_index() const {
    return std::numeric_limits<unsigned>::max();
  }
  constexpr std::vector<unsigned> const &modifiedIndices() const {
    return indices;
  }
};

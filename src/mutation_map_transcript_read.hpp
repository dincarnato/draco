#pragma once

#include <vector>

struct MutationMapTranscriptRead {
  unsigned begin;
  unsigned end;
  std::vector<unsigned> indices;
};

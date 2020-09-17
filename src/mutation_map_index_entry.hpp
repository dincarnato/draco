#pragma once

#include <ios>
#include <string>

struct MutationMapIndexEntry {
  std::string transcriptId;
  std::streampos offset;
};

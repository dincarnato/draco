#pragma once

#include "mutation_map_index_entry.hpp"

#include <string>
#include <string_view>
#include <vector>

struct MutationMapIndex {
  using entries_type = std::vector<MutationMapIndexEntry>;

  MutationMapIndex(std::string_view filename);

  entries_type entries;

private:
  static entries_type readEntries(std::string_view filename);
};

#include "mutation_map_index_impl.hpp"

#pragma once

#include "mutation_map_index.hpp"

#include "binary_stream.hpp"

#include <cassert>

#if __has_include(<filesystem>)
#include <filesystem>
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
#else
#error "Missing filesystem header"
#endif
#include <fstream>

MutationMapIndex::MutationMapIndex(std::string_view filename)
    : entries(readEntries(std::move(filename))) {}

auto MutationMapIndex::readEntries(std::string_view filename) -> entries_type {
#if __has_include(<filesystem>)
  namespace fs = std::filesystem;
#else
  namespace fs = std::experimental::filesystem;
#endif

  std::ifstream indexStream{fs::path(filename)};
  BinaryStream<std::ifstream> binaryStream(indexStream);
  entries_type entries;

  for (;;) {
    std::uint16_t transcriptIdSize;
    binaryStream >> transcriptIdSize;
    if (not indexStream)
      break;

    MutationMapIndexEntry entry;
    entry.transcriptId.resize(transcriptIdSize - 1);
    indexStream.read(entry.transcriptId.data(), transcriptIdSize - 1);

    {
      [[maybe_unused]] int c = indexStream.get();
      assert(c == 0);
    }

    std::uint64_t offset;
    binaryStream >> offset;
    entry.offset = static_cast<std::streamoff>(offset);
    entries.emplace_back(std::move(entry));
  }

  return entries;
}

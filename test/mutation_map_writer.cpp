#include "mutation_map_writer.hpp"

#include "mutation_map.hpp"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <ranges>

int main(int argc, const char *argv[]) {
  assert(argc == 2);
  const char *const filename = argv[1];

  auto output_mm_path = std::filesystem::path(std::tmpnam(nullptr));

  {
    MutationMap mutation_map(filename);
    MutationMapWriter output_mm(output_mm_path);

    for (auto &&transcript : mutation_map) {
      auto output_transcript =
          output_mm.transcript(transcript.getId(), transcript.getSequence());

      for (auto &&read : transcript) {
        output_transcript.add_read(read);
      }
    }
  }

  std::ifstream filename_stream(filename);
  std::ifstream output_mm_stream(output_mm_path);
  assert(
      std::ranges::equal(std::views::istream<std::uint8_t>(filename_stream),
                         std::views::istream<std::uint8_t>(output_mm_stream)));
}

#include "mutation_map_writer.hpp"

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <stdexcept>
#include <string_view>

MutationMapWriter::MutationMapWriter(std::filesystem::path const &path)
    : stream{path, std::ios_base::out | std::ios_base::trunc},
      binary_stream(stream) {}

MutationMapWriter::~MutationMapWriter() {
  if (not stream.is_open()) {
    return;
  }

  stream << "[mmeof]";
  stream.flush();
}

MutationMapTranscriptWriter
MutationMapWriter::transcript(std::string_view id, std::string_view sequence) {
  return MutationMapTranscriptWriter(id, sequence, binary_stream);
}

static std::uint8_t get_base_repr(char base) {
  switch (base) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  case 'N':
    return 4;
  default:
    throw new std::runtime_error("Invalid base");
  }
}

MutationMapTranscriptWriter::MutationMapTranscriptWriter(
    std::string_view id, std::string_view sequence,
    BinaryStream<std::ofstream> &stream)
    : stream(&stream), reads(0) {
  auto id_size = static_cast<std::uint16_t>(id.length() + 1);
  auto sequence_size = static_cast<std::uint16_t>(sequence.length());

  stream << id_size;
  stream.inner() << id << '\0';
  stream << static_cast<std::uint32_t>(sequence_size);

  for (std::uint16_t index = 0; index < sequence_size; index += 2) {
    auto encoded_base =
        static_cast<std::uint8_t>(get_base_repr(sequence[index]) << 4);
    if (index + 1 < sequence_size) {
      encoded_base |= get_base_repr(sequence[index + 1]);
    }
    stream << encoded_base;
  }

  reads_len_pos = stream.tellp();
  stream << static_cast<std::uint32_t>(0);
}

MutationMapTranscriptWriter::~MutationMapTranscriptWriter() {
  auto stream_pos = stream->tellp();
  stream->seekp(reads_len_pos);
  *stream << reads;
  stream->seekp(stream_pos);
}

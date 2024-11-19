#pragma once

#include "binary_stream.hpp"

#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>

struct MutationMapTranscriptWriter;

template <typename T>
concept ReadIndicesLike = requires(T indices) {
  requires std::forward_iterator<
      typename std::remove_reference_t<T>::const_iterator>;
  requires std::is_same_v<typename std::remove_reference_t<T>::value_type,
                          std::uint32_t>;
};

template <typename T>
concept ReadLike = requires(T const &read) {
  { read.original_begin_index() } -> std::same_as<std::uint32_t>;
  { read.original_end_index() } -> std::same_as<std::uint32_t>;
  { read.window_begin_index() } -> std::same_as<std::uint32_t>;
  { read.modifiedIndices() } -> ReadIndicesLike;
};

struct MutationMapWriter {
  explicit MutationMapWriter(std::filesystem::path const &path);
  ~MutationMapWriter();

  MutationMapTranscriptWriter transcript(std::string_view id,
                                         std::string_view sequence);

protected:
  std::ofstream stream;
  BinaryStream<std::ofstream> binary_stream;
};

struct MutationMapTranscriptWriter {
  friend struct MutationMapWriter;

  ~MutationMapTranscriptWriter();

  void add_read(ReadLike auto const &read);

protected:
  explicit MutationMapTranscriptWriter(std::string_view id,
                                       std::string_view sequence,
                                       BinaryStream<std::ofstream> &stream);

  BinaryStream<std::ofstream> *stream;
  std::ofstream::pos_type reads_len_pos;
  std::uint32_t reads;
};

void MutationMapTranscriptWriter::add_read(ReadLike auto const &read) {
  using window_index_t = decltype(read.window_begin_index());

  auto const &modified_indices = read.modifiedIndices();
  *stream << read.original_begin_index() << read.original_end_index() - 1
          << static_cast<uint32_t>(modified_indices.size());

  auto const window_begin_index =
      read.window_begin_index() == std::numeric_limits<window_index_t>::max()
          ? 0
          : read.window_begin_index();
  for (auto index : modified_indices) {
    *stream << index + window_begin_index;
  }

  ++reads;
}

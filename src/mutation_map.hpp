#pragma once

#include "mutation_map_iterator.hpp"
#include "mutation_map_transcript.hpp"
#include "nostd/type_traits.hpp"

#include <array>
#include <deque>
#include <string>

class MutationMapTranscriptHelper;

class MutationMap {
  template <typename> friend class MutationMapIterator;
  friend class MutationMapTranscriptHelper;

public:
  using transcripts_type = std::deque<MutationMapTranscript>;
  using iterator = MutationMapIterator<MutationMap>;
  using const_iterator = MutationMapIterator<const MutationMap>;

  MutationMap() = default;
  MutationMap(std::string_view filename) noexcept(false);

  iterator begin() noexcept;
  iterator end() noexcept;
  const_iterator begin() const noexcept;
  const_iterator end() const noexcept;
  void load() noexcept(false);

  const std::string &getFilename() const noexcept;
  void loadIndexFile() noexcept(false);

private:
  void checkEofMarker() noexcept(false);
  unsigned loadOtherTranscripts(unsigned n) noexcept(false);

  static constexpr std::array<std::uint8_t, 7> eofMarker{
      '\x5b', '\x6d', '\x6d', '\x65', '\x6f', '\x66', '\x5d'};
  std::string filename;
  std::streampos streamEndPos;
  transcripts_type transcripts;
};

#pragma once

#include "mutation_map_transcript_iterator.hpp"
#include "nostd/type_traits.hpp"

#include <fstream>
#include <functional>
#include <string>

class MutationMap;

class MutationMapTranscript {
public:
  using const_iterator = MutationMapTranscriptIterator<std::ifstream>;

  MutationMapTranscript() = default;
  MutationMapTranscript(const MutationMap &mutationMap, const std::string &id,
                        std::streampos offset) noexcept;
  explicit MutationMapTranscript(
      const MutationMap &mutationMap,
      std::streampos startOffset = 0) noexcept(false);

  MutationMapTranscript(const MutationMapTranscript &other) noexcept;
  MutationMapTranscript(MutationMapTranscript &&) = default;
  MutationMapTranscript &
  operator=(const MutationMapTranscript &other) noexcept(false);
  MutationMapTranscript &operator=(MutationMapTranscript &&) = default;

  const std::string &getId() const noexcept;
  const std::string &getSequence() const noexcept;
  unsigned getReadsSize() const noexcept;
  std::streampos getStartOffset() const noexcept;
  std::array<std::vector<std::size_t>, 2> calculateMutationsAndCoverage() const
      noexcept(false);

  void setSequence(const std::string &value) noexcept(false);
  void setSequence(std::string &&value) noexcept(
      std::is_nothrow_move_assignable_v<std::string>);
  void setReads(unsigned value) noexcept;
  bool isLoaded() const noexcept;
  void setEndOffset(const std::streampos &value) const noexcept(false);
  std::streampos getEndOffset() const noexcept;
  const std::string &getFilename() const noexcept;

  void skip() const noexcept(false);
  std::ifstream &getStream() const noexcept(false);

  const_iterator begin() const noexcept(false);
  const_iterator end() const noexcept;

private:
  const MutationMap *mutationMap;
  std::string id;
  std::string sequence;
  unsigned reads;
  std::streampos startOffset;
  mutable std::streampos endOffset;
  mutable std::ifstream iteratorStream;

  void reopenStream() const noexcept(false);
};

#include "mutation_map_transcript_iterator_impl.hpp"

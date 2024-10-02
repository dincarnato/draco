#pragma once

struct MutationMapTranscriptEndTag {};

#include "mutation_map_transcript_read.hpp"

#include <ios>
#include <limits>

class MutationMap;
class MutationMapTranscript;

template <typename Stream> class MutationMapTranscriptIterator {
  friend class MutationMapTranscript;

public:
  using value_type = MutationMapTranscriptRead;
  using reference = value_type const &;
  using pointer = value_type const *;
  using iterator_category = std::input_iterator_tag;
  using difference_type = std::ptrdiff_t;

  using self = MutationMapTranscriptIterator;

  MutationMapTranscriptIterator() = default;
  MutationMapTranscriptIterator(
      const MutationMapTranscript &mutationMapTranscript,
      MutationMapTranscriptEndTag) noexcept;
  explicit MutationMapTranscriptIterator(
      const MutationMapTranscript &mutationMapTranscript,
      unsigned indexOffset = 0) noexcept(false);

  reference operator*() const noexcept;
  pointer operator->() const noexcept;
  self &operator++() noexcept(false);
  self operator++(int) noexcept(false);
  difference_type operator-(const self &other) const noexcept;
  bool operator==(const self &other) const noexcept;
  bool operator!=(const self &other) const noexcept;

private:
  void nextRead() noexcept(false);

  const MutationMapTranscript *mutationMapTranscript = nullptr;
  unsigned index = std::numeric_limits<unsigned>::max();
  value_type currentRead;
};

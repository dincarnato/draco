#pragma once

#include "binary_stream.hpp"
#include "mutation_map_transcript.hpp"
#include "mutation_map_transcript_iterator.hpp"

#include <algorithm>
#include <cassert>

template <typename Stream>
MutationMapTranscriptIterator<Stream>::MutationMapTranscriptIterator(
    const MutationMapTranscript &mutationMapTranscript,
    MutationMapTranscriptEndTag) noexcept
    : mutationMapTranscript(&mutationMapTranscript),
      index(mutationMapTranscript.getReadsSize()) {}

template <typename Stream>
MutationMapTranscriptIterator<Stream>::MutationMapTranscriptIterator(
    const MutationMapTranscript &mutationMapTranscript,
    unsigned indexOffset) noexcept(false)
    : mutationMapTranscript(&mutationMapTranscript), index(indexOffset) {

  if (index < this->mutationMapTranscript->getReadsSize()) {
    for (unsigned index = 0; index < indexOffset; ++index) {
      nextRead();
    }

    nextRead();
  } else {
    this->mutationMapTranscript->setEndOffset(
        this->mutationMapTranscript->getStream().tellg());
  }
}

template <typename Stream>
void MutationMapTranscriptIterator<Stream>::nextRead() noexcept(false) {
  auto &stream = mutationMapTranscript->getStream();
  BinaryStream<Stream> binaryStream(stream);
  {
    std::uint32_t startIndex;
    binaryStream >> startIndex;
    assert(startIndex <= mutationMapTranscript->getSequence().size());
    currentRead.begin = startIndex;
  }

  {
    std::uint32_t endIndex;
    binaryStream >> endIndex;
    assert(endIndex <= mutationMapTranscript->getSequence().size());
    currentRead.end = endIndex + 1;
  }

  {
    std::uint32_t nReads;
    binaryStream >> nReads;
    currentRead.indices.resize(nReads);
    for (auto &readIndex : currentRead.indices) {
      std::uint32_t mutationIndex;
      binaryStream >> mutationIndex;
      readIndex = mutationIndex;
    }

    assert(std::ranges::is_sorted(currentRead.indices));
  }
}

template <typename Stream>
auto MutationMapTranscriptIterator<Stream>::operator*() const noexcept
    -> reference {
  assert(index < mutationMapTranscript->getReadsSize());
  return currentRead;
}

template <typename Stream>
auto MutationMapTranscriptIterator<Stream>::operator->() const noexcept
    -> pointer {
  assert(index < mutationMapTranscript->getReadsSize());
  return &currentRead;
}

template <typename Stream>
auto MutationMapTranscriptIterator<Stream>::operator++() noexcept(false)
    -> self & {
  if (++index < mutationMapTranscript->getReadsSize())
    nextRead();
  else
    mutationMapTranscript->setEndOffset(
        mutationMapTranscript->getStream().tellg());

  return *this;
}

template <typename Stream>
auto MutationMapTranscriptIterator<Stream>::operator++(int) noexcept(false)
    -> self {
  auto copy = *this;
  operator++();
  return copy;
}

template <typename Stream>
auto MutationMapTranscriptIterator<Stream>::operator-(
    const self &other) const noexcept -> difference_type {
  return index - other.index;
}

template <typename Stream>
bool MutationMapTranscriptIterator<Stream>::operator==(
    const self &other) const noexcept {
  return index == other.index;
}

template <typename Stream>
bool MutationMapTranscriptIterator<Stream>::operator!=(
    const self &other) const noexcept {
  return index != other.index;
}

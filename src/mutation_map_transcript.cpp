#include "mutation_map_transcript.hpp"
#include "mutation_map.hpp"

#include <cassert>
#include <iterator>
#include <ranges>

MutationMapTranscript::MutationMapTranscript(
    const MutationMap &mutationMap, std::streampos offset) noexcept(false)
    : mutationMap(&mutationMap), endOffset(-1),
      iteratorStream(mutationMap.getFilename()) {
  assert(offset >= 0);
  iteratorStream.seekg(offset);
  BinaryStream<std::ifstream> binaryStream(iteratorStream);
  {
    std::uint16_t idSizeWithTermination;
    binaryStream >> idSizeWithTermination;
    assert(idSizeWithTermination > 0);

    id.resize(idSizeWithTermination - 1);
    iteratorStream.read(id.data(), idSizeWithTermination - 1);
    {
      [[maybe_unused]] int c = iteratorStream.get();
      assert(c == 0);
    }
  }

  {
    std::uint32_t nBases;
    binaryStream >> nBases;
    assert(nBases > 0);

    sequence.resize(nBases);
    std::vector<char> rawSequence((nBases + 1) / 2);
    iteratorStream.read(rawSequence.data(), (nBases + 1) / 2);

    auto baseIter = std::ranges::begin(sequence);
    for (unsigned baseIndex = 0; baseIndex < (nBases + 1) / 2; ++baseIndex) {
      for (unsigned splitIndex = 0; splitIndex < 2; ++splitIndex, ++baseIter) {
        switch ((rawSequence[baseIndex] >> ((1 - splitIndex) * 4)) & 0x0f) {
        case 0:
          *baseIter = 'A';
          break;
        case 1:
          *baseIter = 'C';
          break;
        case 2:
          *baseIter = 'G';
          break;
        case 3:
          *baseIter = 'T';
          break;
        case 4:
          *baseIter = 'N';
          break;
        }
      }
    }
  }

  {
    std::uint32_t nReads;
    binaryStream >> nReads;
    reads = nReads;
  }

  startOffset = iteratorStream.tellg();
  iteratorStream.close();
}

MutationMapTranscript::MutationMapTranscript(
    const MutationMap &mutationMap, const std::string &id,
    std::streampos startOffset) noexcept
    : mutationMap(&mutationMap), id(id), reads(0), startOffset(startOffset),
      endOffset(-1) {}

MutationMapTranscript::MutationMapTranscript(
    const MutationMapTranscript &other) noexcept
    : mutationMap(other.mutationMap), id(other.id), sequence(other.sequence),
      reads(other.reads), startOffset(other.startOffset),
      endOffset(other.endOffset) {}

MutationMapTranscript &MutationMapTranscript::operator=(
    const MutationMapTranscript &other) noexcept(false) {
  if (other.iteratorStream.is_open()) {
    if (mutationMap != other.mutationMap) {
      if (iteratorStream.is_open())
        iteratorStream.close();

      iteratorStream.open(other.mutationMap->getFilename());
    }
    iteratorStream.seekg(other.startOffset);
  } else if (iteratorStream.is_open())
    iteratorStream.close();

  mutationMap = other.mutationMap;
  id = other.id;
  sequence = other.sequence;
  reads = other.reads;
  startOffset = other.startOffset;
  endOffset = other.endOffset;

  return *this;
}

void MutationMapTranscript::setSequence(const std::string &value) noexcept(
    false) {
  sequence = value;
}

void MutationMapTranscript::setSequence(std::string &&value) noexcept(
    std::is_nothrow_move_assignable_v<std::string>) {
  sequence = std::move(value);
}

void MutationMapTranscript::setReads(unsigned value) noexcept { reads = value; }

const std::string &MutationMapTranscript::getId() const noexcept { return id; }

const std::string &MutationMapTranscript::getSequence() const noexcept {
  return sequence;
}

unsigned MutationMapTranscript::getReadsSize() const noexcept { return reads; }

std::streampos MutationMapTranscript::getStartOffset() const noexcept {
  return startOffset;
}

bool MutationMapTranscript::isLoaded() const noexcept {
  return not sequence.empty();
}

void MutationMapTranscript::setEndOffset(const std::streampos &value) const
    noexcept(false) {
  endOffset = value;

  if (iteratorStream.is_open())
    iteratorStream.close();
}

std::streampos MutationMapTranscript::getEndOffset() const noexcept {
  assert(endOffset != -1);
  return endOffset;
}

const std::string &MutationMapTranscript::getFilename() const noexcept {
  assert(mutationMap);
  return mutationMap->getFilename();
}

auto MutationMapTranscript::begin() const noexcept(false) -> const_iterator {
  return const_iterator(*this);
}

auto MutationMapTranscript::end() const noexcept -> const_iterator {
  return const_iterator(*this, MutationMapTranscriptEndTag{});
}

void MutationMapTranscript::reopenStream() const noexcept(false) {
  if (not iteratorStream.is_open()) {
    iteratorStream.open(mutationMap->getFilename());
    iteratorStream.seekg(startOffset);
  }
}

std::ifstream &MutationMapTranscript::getStream() const noexcept(false) {
  reopenStream();
  return iteratorStream;
}

void MutationMapTranscript::skip() const noexcept(false) {
  if (endOffset != -1) {
    if (iteratorStream.is_open())
      iteratorStream.close();
  } else {
    reopenStream();
    iteratorStream.seekg(startOffset);
    BinaryStream<std::decay_t<decltype(iteratorStream)>> binaryStream(
        iteratorStream);
    for (unsigned readIndex = 0; readIndex < reads; ++readIndex) {
      iteratorStream.seekg(sizeof(std::uint32_t) * 2, std::ios::cur);
      std::uint32_t mutations;
      binaryStream >> mutations;
      iteratorStream.seekg(sizeof(std::uint32_t) * mutations, std::ios::cur);
    }
    endOffset = iteratorStream.tellg();
    iteratorStream.close();
  }
}

std::array<std::vector<std::size_t>, 2>
MutationMapTranscript::calculateMutationsAndCoverage() const noexcept(false) {
  const std::size_t sequenceSize = getSequence().size();
  std::array<std::vector<std::size_t>, 2> out{
      std::vector<std::size_t>(sequenceSize, 0),
      std::vector<std::size_t>(sequenceSize, 0)};

  std::ranges::for_each(*this, [&](auto &&read) {
    for (auto index : read.indices)
      ++out[0][index];

    std::ranges::for_each(out[1] | std::views::drop(read.begin) |
                              std::views::take(read.end - read.begin),
                          [&](auto &index) { ++index; });
  });

  return out;
}

static_assert(std::input_iterator<MutationMapTranscript::const_iterator>);

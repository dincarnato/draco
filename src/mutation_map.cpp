#include "mutation_map.hpp"

#include "mutation_map_index.hpp"
#include "ringmap_data.hpp"

#include <fstream>
#include <range/v3/algorithm.hpp>
#include <range/v3/iterator/insert_iterators.hpp>
#include <string_view>

constexpr std::array<std::uint8_t, 7> MutationMap::eofMarker;

MutationMap::MutationMap(std::string_view filename) : filename(filename) {
  checkEofMarker();
  loadIndexFile();
}

void MutationMap::checkEofMarker() noexcept(false) {
  std::ifstream fileStream(filename, std::ios::binary bitor std::ios::in);
  fileStream.seekg(-static_cast<std::ptrdiff_t>(eofMarker.size()),
                   std::ios::end);
  streamEndPos = fileStream.tellg();

  std::decay_t<decltype(eofMarker)> marker;
  for (auto &value : marker) {
    int stream_char = fileStream.get();
    if (stream_char < 0)
      break;

    value = static_cast<std::uint8_t>(stream_char);
  }
  if (marker != eofMarker)
    throw std::runtime_error("end-of-file magic marker not found");
}

unsigned MutationMap::loadOtherTranscripts(unsigned n) noexcept(false) {
  std::streampos offset;
  if (transcripts.size() > 0) {
    auto &transcript = transcripts.back();
    transcript.skip();
    offset = transcript.getEndOffset();
  } else
    offset = 0;

  unsigned transcriptIndex = 0;
  for (; transcriptIndex < n and offset < streamEndPos; ++transcriptIndex) {
    assert(offset >= 0);
    transcripts.emplace_back(*this, offset);
    if (transcriptIndex < n - 1)
      transcripts.back().skip();
  }
  if (transcripts.size() > 1)
    offset = ranges::prev(ranges::end(transcripts), 2)->getEndOffset();

  return transcriptIndex;
}

auto MutationMap::begin() noexcept -> iterator { return iterator(*this); }

auto MutationMap::end() noexcept -> iterator { return iterator(); }

auto MutationMap::begin() const noexcept -> const_iterator {
  return const_iterator(*this);
}

auto MutationMap::end() const noexcept -> const_iterator {
  return const_iterator();
}

void MutationMap::load() noexcept(false) {
  auto iter = begin();
  while (iter != end())
    ++iter;
}

const std::string &MutationMap::getFilename() const noexcept {
  return filename;
}

void MutationMap::loadIndexFile() noexcept(false) {
  std::string_view filename(this->filename);
  if (filename.substr(std::size(filename) - 3) != ".mm")
    return;

  std::string indexFilename = [&] {
    std::stringstream ss;
    ss << filename << 'i';
    return ss.str();
  }();

  MutationMapIndex mutationMapIndex(indexFilename);
  auto &entries = mutationMapIndex.entries;
  if (entries.empty())
    return;

  ranges::transform(entries, ranges::back_inserter(transcripts),
                    [this](auto &entry) {
                      return MutationMapTranscript(*this, entry.offset);
                    });
}

static_assert(ranges::InputIterator<MutationMap::iterator>);
static_assert(ranges::InputIterator<MutationMap::const_iterator>);

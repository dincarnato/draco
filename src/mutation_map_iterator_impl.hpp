#pragma once

#include "mutation_map_iterator.hpp"

#include <cassert>
#include <stdexcept>

template <typename MMap>
MutationMapIterator<MMap>::MutationMapIterator(MMap& mutationMap,
                                               std::streampos offset) noexcept
    : mutationMap(&mutationMap), offset(offset) {}

template <typename MMap>
auto MutationMapIterator<MMap>::operator*() const noexcept(false) -> reference {
  assert(mutationMap);
  checkOffset();
  return *std::next(std::begin(mutationMap->transcripts), offset);
}

template <typename MMap>
auto MutationMapIterator<MMap>::operator-> () const noexcept(false) -> pointer {
  assert(mutationMap);
  checkOffset();
  return mutationMap->transcripts.data() + offset;
}

template <typename MMap>
auto
MutationMapIterator<MMap>::operator++() noexcept -> self& {
  assert(mutationMap);
  offset += 1;
  return *this;
}

template <typename MMap>
auto
MutationMapIterator<MMap>::operator++(int) noexcept -> self {
  assert(mutationMap);
  self other(*this);
  offset += 1;
  return other;
}

template <typename MMap>
auto
MutationMapIterator<MMap>::operator-(const self& other) const noexcept
    -> difference_type {
  assert(mutationMap);
  return static_cast<difference_type>(offset) - other.offset;
}

template <typename MMap>
template <typename T>
void
MutationMapIterator<MMap>::checkOffset(
    std::enable_if_t<not std::is_const<T>::value>*) const noexcept(false) {
  assert(mutationMap);
  if (static_cast<std::streamoff>(mutationMap->transcripts.size()) <= offset) {
    if (not loadOtherTranscripts())
      throw std::runtime_error("trying to deference a past-over-end iterator");
  }
}

template <typename MMap>
template <typename T>
void
MutationMapIterator<MMap>::checkOffset(
    std::enable_if_t<std::is_const<T>::value>*) const noexcept(false) {
  assert(mutationMap);
  if (static_cast<std::streamoff>(mutationMap->transcripts.size()) <= offset)
    throw std::runtime_error("attempting to deference a not loaded transcript "
                             "on a const mutation map");
}

template <typename MMap>
template <typename T>
std::enable_if_t<not std::is_const<T>::value, bool>
MutationMapIterator<MMap>::loadOtherTranscripts() const noexcept(false) {
  assert(mutationMap);
  auto const transcripts_size =
      static_cast<std::streamoff>(mutationMap->transcripts.size());
  assert(transcripts_size >= offset);
  auto const transcriptsToLoad =
      static_cast<unsigned>(transcripts_size - offset + 1);
  unsigned loadedTranscripts =
      mutationMap->loadOtherTranscripts(transcriptsToLoad);

  return loadedTranscripts == transcriptsToLoad;
}

template <typename MMap>
bool
MutationMapIterator<MMap>::operator==(const self& other) const noexcept(false) {
  if (mutationMap == nullptr) {
    if (other.mutationMap == nullptr)
      return true;

    if (other.offset >=
        static_cast<std::streamoff>(other.mutationMap->transcripts.size()))
      return not other.loadOtherTranscripts();

    return false;
  } else if (other.mutationMap == nullptr) {
    if (offset >= static_cast<std::streamoff>(mutationMap->transcripts.size()))
      return not loadOtherTranscripts();

    return false;
  }

  return offset == other.offset;
}

template <typename MMap>
bool
MutationMapIterator<MMap>::operator!=(const self& other) const noexcept(false) {
  return not operator==(other);
}

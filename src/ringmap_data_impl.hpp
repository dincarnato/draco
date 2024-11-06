#pragma once

#include "ringmap_data.hpp"

template <typename Iter>
RingmapData::RingmapData(const std::string &sequence, unsigned nReads,
                         Iter readsBegin, Iter readsEnd, Args const &args)
    : startIndex(0), endIndex(static_cast<unsigned>(sequence.size())),
      nSetReads(nReads), minimumCoverage(args.minimum_base_coverage()),
      minimumModificationsPerBase(args.minimum_modifications_per_base()),
      minimumModificationsPerRead(args.minimum_modifications_per_read()),
      minimumModificationsPerBaseFraction(
          args.minimum_modifications_per_base_fraction()),
      sequence(sequence), baseCoverages(endIndex, 0), m_data(nReads, endIndex),
      shape(args.shape()) {
  addReads(std::move(readsBegin), std::move(readsEnd));
}

template <typename Iter> void RingmapData::addReads(Iter begin, Iter end) {
#ifndef NDEBUG
  std::size_t readsCount = 0;
#endif
  std::ranges::for_each(begin, end, [&, this](auto &&read) {
    auto beginCurrentBaseCoverages =
        std::ranges::next(std::ranges::begin(baseCoverages), read.begin);
    std::ranges::transform(
        beginCurrentBaseCoverages,
        std::ranges::next(std::ranges::begin(baseCoverages), read.end),
        beginCurrentBaseCoverages, [](unsigned value) { return value + 1; });
    m_data.addRead(read);

#ifndef NDEBUG
    ++readsCount;
#endif
  });

  assert(readsCount == nSetReads);
}

template <typename Iterable>
void RingmapData::keepOnlyIndices(Iterable &&iterable) {
  m_data.keepOnlyIndices(std::forward<Iterable>(iterable));
  nSetReads = m_data.storedReads();
}

template <typename T>
std::enable_if_t<std::is_same<std::decay_t<T>, RingmapData>::value, RingmapData>
RingmapData::from(T &&other) {
  RingmapData ringmap{std::forward<T>(other).sequence,
                      data_type(std::forward<T>(other).m_data),
                      other.startIndex, other.endIndex};
  ringmap.baseCoverages = std::forward<T>(other).baseCoverages;
  ringmap.minimumCoverage = other.minimumCoverage;
  ringmap.minimumModificationsPerBase = other.minimumModificationsPerBase;
  ringmap.minimumModificationsPerRead = other.minimumModificationsPerRead;
  ringmap.minimumModificationsPerBaseFraction =
      other.minimumModificationsPerBaseFraction;

#ifndef NDEBUG
  ringmap.allowedMismatches = other.allowedMismatches;
#endif
  ringmap.shape = other.shape;

  return ringmap;
}

#pragma once

#include <string>

#include "../ringmap_data.hpp"
#include "../ringmap_matrix.hpp"

namespace test {

struct RingmapData {
  static RingmapMatrix const &
  m_data(::RingmapData const &ringmap_data) noexcept {
    return ringmap_data.m_data;
  }

  static RingmapMatrix &m_data(::RingmapData &ringmap_data) noexcept {
    return ringmap_data.m_data;
  }

  static std::string const &
  sequence(::RingmapData const &ringmap_data) noexcept {
    return ringmap_data.sequence;
  }

  static std::string &sequence(::RingmapData &ringmap_data) noexcept {
    return ringmap_data.sequence;
  }

  static unsigned &startIndex(::RingmapData &ringmap_data) noexcept {
    return ringmap_data.startIndex;
  }

  static unsigned startIndex(::RingmapData const &ringmap_data) noexcept {
    return ringmap_data.startIndex;
  }

  static unsigned &endIndex(::RingmapData &ringmap_data) noexcept {
    return ringmap_data.endIndex;
  }

  static unsigned endIndex(::RingmapData const &ringmap_data) noexcept {
    return ringmap_data.endIndex;
  }

  static unsigned &nSetReads(::RingmapData &ringmap_data) noexcept {
    return ringmap_data.nSetReads;
  }

  static unsigned nSetReads(::RingmapData const &ringmap_data) noexcept {
    return ringmap_data.nSetReads;
  }

  static unsigned &minimumCoverage(::RingmapData &ringmap_data) noexcept {
    return ringmap_data.minimumCoverage;
  }

  static unsigned minimumCoverage(::RingmapData const &ringmap_data) noexcept {
    return ringmap_data.minimumCoverage;
  }

  static unsigned &
  minimumModificationsPerBase(::RingmapData &ringmap_data) noexcept {
    return ringmap_data.minimumModificationsPerBase;
  }

  static unsigned
  minimumModificationsPerBase(::RingmapData const &ringmap_data) noexcept {
    return ringmap_data.minimumModificationsPerBase;
  }

  static float minimumModificationsPerBaseFraction(
      ::RingmapData const &ringmap_data) noexcept {
    return ringmap_data.minimumModificationsPerBaseFraction;
  }

  static unsigned
  minimumModificationsPerRead(::RingmapData const &ringmap_data) noexcept {
    return ringmap_data.minimumModificationsPerRead;
  }

  static std::vector<unsigned> const &
  baseCoverages(::RingmapData const &ringmap_data) noexcept {
    return ringmap_data.baseCoverages;
  }
};

} // namespace test

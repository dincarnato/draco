#pragma once

#include "rna_secondary_structure.hpp"

#include <array>
#include <cassert>
#include <limits>
#include <vector>

struct Helix {
  using iter_type = std::size_t;
  using offset_type = std::make_signed_t<iter_type>;

  struct PairsSpan {
    iter_type begin;
    iter_type end;
    std::size_t pairingSize;
  };

  Helix() = default;
  Helix(iter_type begin, iter_type end,
        const RnaSecondaryStructure& structure) noexcept
      : _begin(std::move(begin)), _end(std::move(end)) {
    setPairs(structure);
  }

  explicit Helix(iter_type begin) noexcept
      : _begin(std::move(begin))
#ifndef NDEBUG
        ,
        _end(std::numeric_limits<iter_type>::max())
#endif
  {
  }

  bool
  isCompatible(const Helix& other) const noexcept {
    assert(_end != std::numeric_limits<iter_type>::max());
    assert(other._end != std::numeric_limits<iter_type>::max());
    assert(not _pairs.empty());
    assert(not other._pairs.empty());

    if (_begin >= other._end or _end <= other._begin)
      return true;

    if ((_begin >= other._begin and _begin < other._end and
         _end >= other._end) or
        (other._begin >= _begin and other._begin < _end and other._end >= _end))
      return false;

    for (const auto& pair : _pairs) {
      for (const auto& otherPair : other._pairs) {
        if (not(other._begin + otherPair.begin >= _begin + pair.end or
                (_begin + pair.begin >=
                     other._begin + otherPair.begin + otherPair.pairingSize and
                 _begin + pair.end <=
                     other._begin + otherPair.end - otherPair.pairingSize)) and
            not(_begin + pair.begin >= other._begin + otherPair.end or
                (other._begin + otherPair.begin >=
                     _begin + pair.begin + pair.pairingSize and
                 other._begin + otherPair.end <=
                     _begin + pair.end - pair.pairingSize)))
          return false;
      }
    }

    return true;
  }

  void
  setPairs(const RnaSecondaryStructure& structure) {
    const auto beginIter =
        std::next(std::begin(structure), static_cast<offset_type>(_begin));
    auto fwIter = beginIter;
    auto bwIter =
        std::next(std::begin(structure), static_cast<offset_type>(_end)) - 1;

    _pairs.clear();

    PairsSpan span{0, _end - _begin, 0};
    while (fwIter < bwIter) {
      assert(fwIter < std::end(structure));
      assert(bwIter < std::end(structure));
      assert(fwIter >= beginIter);
      assert(bwIter >= beginIter);
      if (*fwIter == BaseSecondaryStructure::double_strand_close or
          *bwIter == BaseSecondaryStructure::double_strand_open)
        break;

      if (*fwIter == BaseSecondaryStructure::double_strand_open and
          *bwIter == BaseSecondaryStructure::double_strand_close) {
        ++span.pairingSize;
        ++fwIter;
        --bwIter;
      } else {
        if (span.pairingSize > 0) {
          _pairs.emplace_back(span);
          span.pairingSize = 0;
        }

        while (fwIter != std::end(structure) and
               *fwIter == BaseSecondaryStructure::single_strand)
          ++fwIter;
        span.begin = static_cast<iter_type>(std::distance(beginIter, fwIter));

        while (bwIter != beginIter and
               *bwIter == BaseSecondaryStructure::single_strand)
          --bwIter;
        span.end = static_cast<iter_type>(std::distance(beginIter, bwIter + 1));

        assert(span.begin < _end - _begin);
        assert(span.end < _end - _begin);
      }
    }

    if (span.pairingSize > 0)
      _pairs.emplace_back(std::move(span));
  }

  iter_type
  begin() const noexcept {
    return _begin;
  }

  iter_type
  end() const noexcept {
    return _end;
  }

  void
  setEnd(iter_type end, const RnaSecondaryStructure& structure) {
    _end = end;
    setPairs(structure);
  }

  bool
  operator<(const Helix& other) const noexcept {
    if (_begin < other._begin)
      return true;
    else if (_begin > other._begin)
      return false;

    if (_end < other._end)
      return true;
    else if (_end > other._end)
      return false;

    if (_pairs.size() < other._pairs.size())
      return true;
    else if (_pairs.size() > other._pairs.size())
      return false;

    for (auto pairsIter = std::begin(_pairs),
              otherPairsIter = std::begin(other._pairs);
         pairsIter != std::end(_pairs); ++pairsIter, ++otherPairsIter) {
      const auto& pair = *pairsIter;
      const auto& otherPair = *otherPairsIter;

      if (pair.begin < otherPair.begin)
        return true;
      if (pair.begin > otherPair.begin)
        return false;

      if (pair.pairingSize < otherPair.pairingSize)
        return true;
      if (pair.pairingSize > otherPair.pairingSize)
        return false;

      if (pair.end < otherPair.end)
        return true;
      else if (pair.end > otherPair.end)
        return false;
    }

    return false;
  }

  const std::vector<PairsSpan>&
  pairs() const noexcept {
    return _pairs;
  }

  void
  shift(offset_type offset) {
    assert(offset < 0 ? _begin >= static_cast<iter_type>(-offset) : true);
    _begin = static_cast<iter_type>(static_cast<offset_type>(_begin) + offset);
    _end += static_cast<iter_type>(static_cast<offset_type>(_end) + offset);
  }

private:
  iter_type _begin;
  iter_type _end;
  std::vector<PairsSpan> _pairs;
};

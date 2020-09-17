#include "paired_rna_secondary_structure.hpp"

#include <algorithm>
#include <stack>
#include <stdexcept>

auto
PairedRnaSecondaryStructure::createPairs() const -> pairs_type {
  pairs_type pairs(structure.size(),
                   PairedBase{std::numeric_limits<std::size_t>::max()});

  std::stack<std::size_t> boundBases;

  for (std::size_t baseIndex = 0; baseIndex < structure.size(); ++baseIndex) {
    switch (structure[baseIndex]) {
    case BaseSecondaryStructure::single_strand:
      continue;
    case BaseSecondaryStructure::double_strand_open:
      boundBases.emplace(baseIndex);
      break;
    case BaseSecondaryStructure::double_strand_close:
      if (boundBases.empty())
        throw std::runtime_error("unbalanced structure");

      std::size_t otherBaseIndex = boundBases.top();
      pairs[baseIndex].index = otherBaseIndex;
      pairs[otherBaseIndex].index = baseIndex;

      boundBases.pop();
      break;
    }
  }

  if (not boundBases.empty())
    throw std::runtime_error("unbalanced structure");

  return pairs;
}

std::size_t
PairedRnaSecondaryStructure::countEqualPairs(
    const PairedRnaSecondaryStructure& other) const {
  if (structure.size() != other.structure.size())
    throw std::runtime_error(
        "cannot compare structures with different lengths");

  std::size_t count = 0;

  const auto thisEnd = std::end(pairs);
  auto thisIter = std::begin(pairs);
  auto otherIter = std::begin(other.pairs);
  for (unsigned baseIndex = 0; thisIter != thisEnd;
       ++thisIter, ++otherIter, ++baseIndex) {
    const auto& paired = *thisIter;
    if (not paired.isPaired() or baseIndex > paired())
      continue;

    const auto& otherPaired = *otherIter;
    if (paired() == otherPaired())
      ++count;
  }

  return count;
}

auto
PairedRnaSecondaryStructure::calculateStatistics(
    const PairedRnaSecondaryStructure& predicted,
    const PairedRnaSecondaryStructure& expected) -> Statistics {
  const auto countPaired = [](const auto& structure) {
    return structure.structure.size() -
           static_cast<std::size_t>(std::count(
               std::begin(structure.structure), std::end(structure.structure),
               BaseSecondaryStructure::single_strand));
  };

  const double ppv = static_cast<double>(predicted.countEqualPairs(expected)) *
                     2. / static_cast<double>(countPaired(predicted));
  const double sensitivity =
      static_cast<double>(expected.countEqualPairs(predicted)) * 2. /
      static_cast<double>(countPaired(expected));
  return Statistics{ppv, sensitivity};
}

std::size_t
distance(const PairedRnaSecondaryStructure& structureA,
         const PairedRnaSecondaryStructure& structureB) {
  const auto size = structureA.structure.size();
  assert(structureA.pairs.size() == size);
  assert(structureB.pairs.size() == size);
  assert(structureB.structure.size() == size);
  std::size_t distance = size;

  {
    auto pairsAIter = std::begin(structureA.pairs);
    auto pairsBIter = std::begin(structureB.pairs);
    const auto pairsAEnd = std::end(structureA.pairs);
    for (; pairsAIter != pairsAEnd; ++pairsAIter, ++pairsBIter) {
      if (pairsAIter->index == pairsBIter->index)
        --distance;
    }
  }

  return distance;
}

const RnaSecondaryStructure&
PairedRnaSecondaryStructure::raw() const {
  return structure;
}

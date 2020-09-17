#pragma once

#include "paired_rna_secondary_structure.hpp"

#include <cassert>
#include <limits>

inline bool
PairedBase::isPaired() const {
  return index != std::numeric_limits<std::size_t>::max();
}

inline std::size_t
PairedBase::operator()() const {
  return index;
}

inline bool
PairedBase::operator<(const PairedBase& other) const {
  return index < other.index;
}

inline PairedRnaSecondaryStructure::PairedRnaSecondaryStructure(
    const RnaSecondaryStructure& structure)
    : structure(structure), pairs(createPairs()) {}

inline PairedRnaSecondaryStructure::PairedRnaSecondaryStructure(
    RnaSecondaryStructure&& structure)
    : structure(std::move(structure)), pairs(createPairs()) {}

inline bool
PairedRnaSecondaryStructure::isPaired(std::size_t index) const {
  assert(index < pairs.size());
  return pairs[index].isPaired();
}

inline const PairedBase&
PairedRnaSecondaryStructure::paired(std::size_t index) const {
  assert(index < pairs.size());
  return pairs[index];
}

inline auto
PairedRnaSecondaryStructure::begin() const -> const_iterator {
  return std::begin(pairs);
}

inline auto
PairedRnaSecondaryStructure::end() const -> const_iterator {
  return std::end(pairs);
}

inline const PairedBase& PairedRnaSecondaryStructure::
operator[](std::size_t index) const {
  return pairs[index];
}

inline const RnaSecondaryStructure&
PairedRnaSecondaryStructure::getSecondaryStructure() const {
  return structure;
}

inline auto
PairedRnaSecondaryStructure::size() const -> size_type {
  return structure.size();
}

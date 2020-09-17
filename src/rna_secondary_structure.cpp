#include "rna_secondary_structure.hpp"
#include "helix.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <stack>
#include <stdexcept>

#include <iostream>

RnaSecondaryStructure::RnaSecondaryStructure(const std::string& rawStructure)
    : base_type(), rawStructure(rawStructure) {
  base_type::reserve(rawStructure.size());
  for (char baseRawStructure : rawStructure) {
    if (baseRawStructure == '.' or baseRawStructure == ',')
      base_type::push_back(BaseSecondaryStructure::single_strand);
    else if (baseRawStructure == '(')
      base_type::push_back(BaseSecondaryStructure::double_strand_open);
    else if (baseRawStructure == ')')
      base_type::push_back(BaseSecondaryStructure::double_strand_close);
    else
      throw std::runtime_error("invalid dot bracket");
  }
}

RnaSecondaryStructure::RnaSecondaryStructure(const char* const rawStructure)
    : base_type(), rawStructure(rawStructure) {
  base_type::reserve(std::strlen(rawStructure));
  for (const char* baseRawStructure = rawStructure; *baseRawStructure != '\0';
       ++baseRawStructure) {
    if (*baseRawStructure == '.' or *baseRawStructure == ',')
      base_type::push_back(BaseSecondaryStructure::single_strand);
    else if (*baseRawStructure == '(')
      base_type::push_back(BaseSecondaryStructure::double_strand_open);
    else if (*baseRawStructure == ')')
      base_type::push_back(BaseSecondaryStructure::double_strand_close);
    else
      throw std::runtime_error("invalid dot bracket");
  }
}

RnaSecondaryStructure::RnaSecondaryStructure(std::size_t length)
    : base_type(length), rawStructure(length, '.') {}

RnaSecondaryStructure::RnaSecondaryStructure(
    const RnaSecondaryStructure& other) noexcept
    : base_type(other), rawStructure(other.rawStructure) {}

RnaSecondaryStructure::RnaSecondaryStructure(
    RnaSecondaryStructure&& other) noexcept
    : base_type(std::move(other)), rawStructure(std::move(other.rawStructure)) {
}

unsigned
RnaSecondaryStructure::distance(const RnaSecondaryStructure& other) const {
  auto thisIter = std::begin(*this);
  auto otherIter = std::begin(other);

  unsigned distance = 0;
  for (; thisIter < std::end(*this); ++thisIter, ++otherIter) {
    if (*thisIter != *otherIter)
      ++distance;
  }

  return distance;
}

const std::string&
RnaSecondaryStructure::raw() const {
  return rawStructure;
}

RnaSecondaryStructure&
RnaSecondaryStructure::operator=(const RnaSecondaryStructure& other) noexcept {
  base_type::operator=(other);
  rawStructure = other.rawStructure;

  return *this;
}

RnaSecondaryStructure&
RnaSecondaryStructure::operator=(RnaSecondaryStructure&& other) noexcept {
  base_type::operator=(std::move(other));
  rawStructure = std::move(other.rawStructure);

  return *this;
}

RnaSecondaryStructure&
RnaSecondaryStructure::operator+=(const RnaSecondaryStructure& other) {
  std::size_t oldSize = size();
  resize(oldSize + other.size());
  std::copy(std::begin(other), std::end(other),
            std::next(std::begin(*this), static_cast<std::ptrdiff_t>(oldSize)));
  rawStructure += other.rawStructure;

  return *this;
}

RnaSecondaryStructure&
RnaSecondaryStructure::operator+=(RnaSecondaryStructure&& other) {
  std::size_t oldSize = size();
  resize(oldSize + other.size());
  std::move(std::begin(other), std::end(other),
            std::next(std::begin(*this), static_cast<std::ptrdiff_t>(oldSize)));
  rawStructure += std::move(other.rawStructure);

  return *this;
}

RnaSecondaryStructure
RnaSecondaryStructure::operator+(const RnaSecondaryStructure& other) const {
  RnaSecondaryStructure out(*this);
  return out += other;
}

RnaSecondaryStructure
RnaSecondaryStructure::operator+(RnaSecondaryStructure&& other) const {
  RnaSecondaryStructure out(*this);
  return out += std::move(other);
}

std::string
RnaSecondaryStructure::str() const {
  auto strandness = std::begin(*this);
  std::string rawStructure;
  rawStructure.resize(size());

  auto rawStrandness = std::begin(rawStructure);
  for (; strandness < std::end(*this); ++strandness, ++rawStrandness) {
    if (*strandness == BaseSecondaryStructure::single_strand)
      *rawStrandness = '.';
    else if (*strandness == BaseSecondaryStructure::double_strand_open)
      *rawStrandness = '(';
    else if (*strandness == BaseSecondaryStructure::double_strand_close)
      *rawStrandness = ')';
  }

  return rawStructure;
}

void
RnaSecondaryStructure::regenRaw() {
  rawStructure = str();
}

bool
RnaSecondaryStructure::isBalanced() const noexcept {
  return std::count(std::begin(*this), std::end(*this),
                    BaseSecondaryStructure::double_strand_open) ==
         std::count(std::begin(*this), std::end(*this),
                    BaseSecondaryStructure::double_strand_close);
}

std::vector<Helix>
RnaSecondaryStructure::getHelices(bool keepLonelyPairs) const {
  std::vector<Helix> helices;
  auto structure = *this;

  for (;;) {
    assert(structure.isBalanced());
    auto closingIter = std::find(std::begin(structure), std::end(structure),
                                 BaseSecondaryStructure::double_strand_close);
    if (closingIter == std::end(structure))
      break;

    auto openingIter =
        std::find(std::make_reverse_iterator(closingIter), std::rend(structure),
                  BaseSecondaryStructure::double_strand_open);
    assert(openingIter != std::rend(structure));

    for (; closingIter != std::end(structure) and
           openingIter != std::rend(structure);) {
      const auto closing = *closingIter;
      const auto opening = *openingIter;

      if (closing == BaseSecondaryStructure::double_strand_open or
          opening == BaseSecondaryStructure::double_strand_close or
          (closing == BaseSecondaryStructure::single_strand and
           opening == closing))
        break;

      if (closing == BaseSecondaryStructure::double_strand_close and
          opening == BaseSecondaryStructure::double_strand_open) {
        ++closingIter;
        ++openingIter;
      } else {
        while (openingIter != std::rend(structure) and
               *openingIter == BaseSecondaryStructure::single_strand)
          ++openingIter;

        while (closingIter != std::end(structure) and
               *closingIter == BaseSecondaryStructure::single_strand)
          ++closingIter;
      }
    }

    auto startingIter = openingIter.base();

    if (startingIter == std::begin(structure) or
        *startingIter == BaseSecondaryStructure::single_strand)
      startingIter = std::find(startingIter, closingIter,
                               BaseSecondaryStructure::double_strand_open);
    assert(startingIter >= std::begin(structure) and
           startingIter <= std::end(structure));

    if (*std::prev(closingIter) == BaseSecondaryStructure::single_strand)
      closingIter =
          std::find(std::make_reverse_iterator(closingIter), openingIter,
                    BaseSecondaryStructure::double_strand_close)
              .base();
    assert(closingIter >= std::begin(structure) and
           closingIter <= std::end(structure));

    // std::cout << "staringIter=" << std::distance(std::begin(structure),
    // startingIter) << " closingIter=" << std::distance(std::begin(structure),
    // closingIter) << std::endl;
    assert(std::count(startingIter, closingIter,
                      BaseSecondaryStructure::double_strand_open) ==
           std::count(startingIter, closingIter,
                      BaseSecondaryStructure::double_strand_close));

    assert(*startingIter == BaseSecondaryStructure::double_strand_open);
    assert(*std::prev(closingIter) ==
           BaseSecondaryStructure::double_strand_close);

    if (keepLonelyPairs)
      helices.emplace_back(std::distance(std::begin(structure), startingIter),
                           std::distance(std::begin(structure), closingIter),
                           structure);
    else {
      Helix helix(static_cast<std::size_t>(
                      std::distance(std::begin(structure), startingIter)),
                  static_cast<std::size_t>(
                      std::distance(std::begin(structure), closingIter)),
                  structure);
      if (helix.pairs().size() > 1 or helix.pairs()[0].pairingSize > 1)
        helices.push_back(std::move(helix));
    }

    std::fill(startingIter, closingIter, BaseSecondaryStructure::single_strand);
  }

  assert(std::all_of(
      std::begin(structure), std::end(structure), [](const auto strandness) {
        return strandness == BaseSecondaryStructure::single_strand;
      }));

  return helices;
}

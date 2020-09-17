#pragma once

#include "rna_secondary_structure.hpp"

struct PairedBase {
  std::size_t index;

  bool isPaired() const;
  std::size_t operator()() const;

  bool operator<(const PairedBase& other) const;
};

struct PairedRnaSecondaryStructure {
  using pairs_type = std::vector<PairedBase>;
  using const_iterator = typename pairs_type::const_iterator;
  using size_type = typename RnaSecondaryStructure::size_type;

  struct Statistics {
    double ppv;
    double sensitivity;
  };

  PairedRnaSecondaryStructure() = default;
  explicit PairedRnaSecondaryStructure(const RnaSecondaryStructure& structure);
  explicit PairedRnaSecondaryStructure(RnaSecondaryStructure&& structure);

  bool isPaired(std::size_t index) const;
  const PairedBase& paired(std::size_t index) const;

  const_iterator begin() const;
  const_iterator end() const;
  const PairedBase& operator[](std::size_t index) const;

  const RnaSecondaryStructure& getSecondaryStructure() const;
  std::size_t countEqualPairs(const PairedRnaSecondaryStructure& other) const;

  static Statistics
  calculateStatistics(const PairedRnaSecondaryStructure& predicted,
                      const PairedRnaSecondaryStructure& exptected);

  friend std::size_t distance(const PairedRnaSecondaryStructure&,
                              const PairedRnaSecondaryStructure&);
  const RnaSecondaryStructure& raw() const;
  size_type size() const;

private:
  RnaSecondaryStructure structure;
  pairs_type pairs;

  pairs_type createPairs() const;
};

std::size_t distance(const PairedRnaSecondaryStructure&,
                     const PairedRnaSecondaryStructure&);

#include "paired_rna_secondary_structure_impl.hpp"

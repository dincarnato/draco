#include "ringmap_data.hpp"
#include "mutation_map.hpp"
#include "mutation_map_transcript.hpp"
#include "ptba.hpp"
#include "results/window.hpp"
#include "rna_secondary_structure.hpp"
#include "spectral_partitioner.hpp"
#include "tokenizer_iterator.hpp"

#include <range/v3/algorithm.hpp>
#include <range/v3/iterator/insert_iterators.hpp>
#include <range/v3/numeric.hpp>

#include <array>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <regex>
#include <sstream>
#include <unordered_map>

RingmapData::RingmapData(const std::string &filename,
                         const std::string &sequence, Args const &args,
                         bool keepFragments)
    : minimumCoverage(args.minimum_base_coverage()),
      minimumModificationsPerBase(args.minimum_modifications_per_base()),
      minimumModificationsPerRead(args.minimum_modifications_per_read()),
      minimumModificationsPerBaseFraction(
          args.minimum_modifications_per_base_fraction()),
      sequence(sequence), shape(args.shape()) {
  startIndex = std::numeric_limits<unsigned>::max();
  endIndex = std::numeric_limits<unsigned>::lowest();

  std::ifstream stream(filename);
  for (std::string line; std::getline(stream, line);) {
    unsigned currentStartIndex, currentEndIndex;
    std::stringstream lineBuffer(line);
    lineBuffer >> currentStartIndex >> currentEndIndex;

    if (currentStartIndex < startIndex)
      startIndex = currentStartIndex;

    if (currentEndIndex > endIndex)
      endIndex = currentEndIndex;

    ++nSetReads;
  }

  --startIndex;
  m_data = data_type(nSetReads, endIndex - startIndex);
  stream.clear();
  stream.seekg(0);
  baseCoverages.resize(sequence.size(), 0);

  unsigned readIndex = 0;
  for (std::string line; std::getline(stream, line);) {
    auto tokenizer = makeTokenizerIterator(line, '\t');

    unsigned currentStartIndex = tokenizer++.tou() - 1;
    unsigned currentEndIndex = tokenizer++.tou();
    if (not keepFragments) {
      if (currentStartIndex != startIndex or currentEndIndex != endIndex)
        continue;
    }

    {
      auto beginCurrentBaseCoverages =
          ranges::next(ranges::begin(baseCoverages), currentStartIndex);
      ranges::transform(
          beginCurrentBaseCoverages,
          ranges::next(ranges::begin(baseCoverages), currentEndIndex),
          beginCurrentBaseCoverages, [](unsigned value) { return value + 1; });
    }

    const auto rawHits = *tokenizer;
    auto &&readHits = m_data.row(readIndex);

    unsigned hitIndex = 0;
    for (char hit : rawHits) {
      if (hit == '1')
        readHits[currentStartIndex + hitIndex] = 1;

      ++hitIndex;
    }

    ++readIndex;
  }

  if (readIndex != nSetReads)
    m_data.remove_rows(readIndex, nSetReads);
  nSetReads = readIndex;
}

RingmapData::RingmapData(const std::string &sequence, data_type &&dataMatrix,
                         unsigned startIndex, unsigned endIndex,
                         Args const &args)
    : startIndex(startIndex), endIndex(endIndex),
      nSetReads(dataMatrix.storedReads()),
      minimumCoverage(args.minimum_base_coverage()),
      minimumModificationsPerBase(args.minimum_modifications_per_base()),
      minimumModificationsPerRead(args.minimum_modifications_per_read()),
      minimumModificationsPerBaseFraction(
          args.minimum_modifications_per_base_fraction()),
      sequence(sequence), baseCoverages(sequence.size(), nSetReads),
      m_data(std::move(dataMatrix)), shape(args.shape()) {}

RingmapData::RingmapData(const std::string &sequence, data_type &&dataMatrix,
                         unsigned startIndex, unsigned endIndex)
    : startIndex(startIndex), endIndex(endIndex),
      nSetReads(dataMatrix.storedReads()), sequence(sequence),
      baseCoverages(sequence.size(), nSetReads), m_data(std::move(dataMatrix)) {
}

RingmapData::RingmapData(const MutationMapTranscript &transcript,
                         Args const &args)
    : startIndex(0),
      endIndex(static_cast<unsigned>(transcript.getSequence().size())),
      nSetReads(transcript.getReadsSize()),
      minimumCoverage(args.minimum_base_coverage()),
      minimumModificationsPerBase(args.minimum_modifications_per_base()),
      minimumModificationsPerRead(args.minimum_modifications_per_read()),
      minimumModificationsPerBaseFraction(
          args.minimum_modifications_per_base_fraction()),
      sequence(transcript.getSequence()), baseCoverages(endIndex),
      m_data(nSetReads, endIndex), shape(args.shape()) {
  addReads(ranges::begin(transcript), ranges::end(transcript));
}

RingmapData::RingmapData(const RingmapData &other,
                         const std::vector<unsigned> &subsetIndices)
    :
#ifndef NDEBUG
      allowedMismatches(other.allowedMismatches),
#endif
      basesFiltered(other.basesFiltered), startIndex(other.startIndex),
      endIndex(other.endIndex),
      nSetReads(static_cast<unsigned>(subsetIndices.size())),
      modificationsFilter(other.modificationsFilter),
      minimumCoverage(other.minimumCoverage),
      minimumModificationsPerBase(other.minimumModificationsPerBase),
      minimumModificationsPerRead(other.minimumModificationsPerRead),
      minimumModificationsPerBaseFraction(
          other.minimumModificationsPerBaseFraction),
      sequence(other.sequence), baseCoverages(other.baseCoverages),
      oldColsToNew(other.oldColsToNew),
      m_data(static_cast<unsigned>(subsetIndices.size()),
             other.m_data.cols_size()),
      shape(other.shape) {
  m_data.resize(static_cast<unsigned>(subsetIndices.size()));
  unsigned row = 0;
  for (unsigned index : subsetIndices)
    m_data.row(row++) = other.m_data.row(index);
}

void RingmapData::filterBases() {
  assert(basesFiltered or
         ranges::all_of(baseCoverages,
                        [nReads = m_data.rows_size()](unsigned coverage) {
                          return coverage <= nReads;
                        }));

  std::vector<std::size_t> allowedBases;
  if (basesMask.empty() or basesFiltered) {
    allowedBases.resize(m_data.cols_size());
    ranges::iota(allowedBases, static_cast<std::size_t>(0));
  } else {
    allowedBases.resize(static_cast<std::size_t>(ranges::count(basesMask, 1.)));
    auto allowedBaseIter = ranges::begin(allowedBases);
    std::size_t baseIndex = 0;

    for (bool allowed : basesMask) {
      if (allowed)
        *allowedBaseIter++ = baseIndex;
      ++baseIndex;
    }

    assert(allowedBaseIter == ranges::end(allowedBases));
  }

  std::vector<unsigned> usedCols;
  usedCols.reserve(allowedBases.size());

  {
    data_type filtered(m_data.storedReads(),
                       static_cast<unsigned>(allowedBases.size()));
    unsigned new_n_cols = 0;

    decltype(oldColsToNew) newOldColsToNew;

    for (std::size_t allowedIndex = 0; allowedIndex < allowedBases.size();
         ++allowedIndex) {

      std::size_t colIndex = allowedBases[allowedIndex];
      std::size_t baseIndex = [&, this] {
        if (basesFiltered) {
          auto oldColToNewIter =
              ranges::find_if(oldColsToNew, [&](const auto oldAndNew) {
                return oldAndNew.second == colIndex;
              });

          assert(oldColToNewIter != ranges::end(oldColsToNew));
          return static_cast<std::size_t>(oldColToNewIter->first);
        } else
          return colIndex;
      }();

      char base = static_cast<char>(
          std::toupper(static_cast<int>(sequence[baseIndex])));
      assert(not basesFiltered or shape or base == 'C' or base == 'A');

      if (shape or base == 'C' or base == 'A') {
        auto takeCol = [&] {
          usedCols.emplace_back(colIndex);
          newOldColsToNew.emplace(baseIndex, new_n_cols);
          filtered.col(new_n_cols++) =
              m_data.col(static_cast<unsigned>(colIndex));
        };

        auto const modificationsOnCol = static_cast<double>(
            m_data.col(static_cast<unsigned>(colIndex)).sum());
        if (not basesFiltered) {
          if (baseCoverages[colIndex] >= minimumCoverage and
              modificationsOnCol / baseCoverages[colIndex] >
                  minimumModificationsPerBaseFraction and
              modificationsOnCol > minimumModificationsPerBase)
            takeCol();
        } else {
          if (modificationsOnCol > minimumModificationsPerBase)
            takeCol();
        }
      }
    }

    if (oldColsToNew.empty() or new_n_cols != filtered.cols_size()) {
      oldColsToNew = std::move(newOldColsToNew);
      filtered.remove_cols(new_n_cols, filtered.cols_size());
      m_data = std::move(filtered);
      baseCoverages.erase(ranges::transform(usedCols,
                                            ranges::begin(baseCoverages),
                                            [this](std::size_t index) {
                                              return baseCoverages[index];
                                            })
                              .out,
                          ranges::end(baseCoverages));
    }
  }

  basesFiltered = true;
}

void RingmapData::filterReads() {
  if (modificationsFilter >= minimumModificationsPerRead)
    return;

  if (readsMap.empty()) {
    readsMap.resize(m_data.storedReads());
    ranges::iota(readsMap, 0u);
  }

  decltype(readsMap) newReadsMap;
  newReadsMap.reserve(readsMap.size());

  data_type filtered(m_data.storedReads(), m_data.cols_size());
  unsigned new_n_rows = 0;

  for (unsigned row = 0; row < m_data.rows_size(); ++row) {
    if (m_data.row(row).sum() >= minimumModificationsPerRead) {
      newReadsMap.emplace_back(readsMap[row]);
      filtered.row(new_n_rows++) = m_data.row(row);
    }
  }

  filtered.remove_rows(new_n_rows, filtered.rows_size());
  m_data = std::move(filtered);
  readsMap = std::move(newReadsMap);
  nSetReads = new_n_rows;

  modificationsFilter = minimumModificationsPerRead;
}

void RingmapData::filter() {
  filterBases();
  filterReads();
}

void RingmapData::perturb() {
  std::mt19937 randomGenerator(std::random_device{}());
  for (auto &&col : m_data.cols())
    col.shuffle(randomGenerator);
}

auto RingmapData::data() const -> const data_type & { return m_data; }

void RingmapData::removeHighValuesOnAdjacency(arma::mat &adjacency,
                                              double maxFraction) {
  for (unsigned col = 0; col < adjacency.n_cols; ++col) {
    double diagValue = adjacency(col, col);
    for (unsigned row = col + 1; row < adjacency.n_rows; ++row) {
      double &value = adjacency(row, col);
      if (value > maxFraction * diagValue) {
        value = 0.;
        adjacency(col, row) = 0.;
      }
    }
  }
}

bool RingmapData::isFiltered() const {
  return basesFiltered and modificationsFilter > 0;
}

auto RingmapData::getNonFilteredToFilteredMap() const -> const
    decltype(oldColsToNew) & {
  assert(basesFiltered);
  return oldColsToNew;
}

auto RingmapData::getFilteredToNonFilteredMap() const
    -> decltype(oldColsToNew) {
  assert(basesFiltered);
  decltype(oldColsToNew) newColsToOld;
  ranges::transform(oldColsToNew,
                    ranges::inserter(newColsToOld, ranges::begin(newColsToOld)),
                    [](const auto &oldAndNew) {
                      return std::make_pair(oldAndNew.second, oldAndNew.first);
                    });
  return newColsToOld;
}

auto RingmapData::getReadsMap() const -> const decltype(readsMap) & {
  assert(modificationsFilter > 0);
  return readsMap;
}

void RingmapData::shed(unsigned begin, unsigned end) {
  assert(not basesFiltered);
  assert(modificationsFilter > 0);

  end = std::min(end, m_data.cols_size());
  sequence =
      sequence.substr(0, begin) + sequence.substr(end, sequence.size() - end);

  m_data.remove_cols(begin, end);
}

void RingmapData::remove(unsigned begin, unsigned end) {
  end = std::min(end, m_data.cols_size());

  /* We assume to change the baseCoverage proportionally */
  ranges::transform(
      baseCoverages, ranges::begin(baseCoverages),
      [oldSize = m_data.rows_size(),
       newSize = m_data.rows_size() + begin - end](auto coverage) {
        return static_cast<unsigned>(static_cast<double>(coverage) / oldSize *
                                     newSize);
      });

  if (readsMap.size() > 0)
    readsMap.erase(ranges::next(ranges::begin(readsMap), begin),
                   ranges::next(ranges::begin(readsMap), end));
  m_data.remove_rows(begin, end);
  nSetReads -= end - begin;

  assert(nSetReads == m_data.rows_size());
}

void RingmapData::dump(const std::string &filename) const {
  std::ofstream outputFile(filename);
  for (auto &&rowIter = m_data.rows().begin(); rowIter < m_data.rows().end();
       ++rowIter) {
    auto &&row = *rowIter;
    outputFile << "1\t" << sequence.size() << '\t';
    for (auto &&valueIter = row.begin(); valueIter < row.end(); ++valueIter)
      std::cout << *valueIter << ' ';
    std::cout << '\n';
  }
}

RingmapData &RingmapData::operator+=(const RingmapData &other) {
  assert(startIndex == other.startIndex);
  assert(endIndex == other.endIndex);

#ifndef NDEBUG
  if (not allowedMismatches)
    assert(sequence == other.sequence);
#endif

  assert(not basesFiltered);
  assert(not other.basesFiltered);
  if (modificationsFilter != 0) {
    modificationsFilter = 0;
    readsMap.clear();
  }

  m_data.append(other.m_data);
  nSetReads += other.nSetReads;
  ranges::transform(other.baseCoverages, baseCoverages,
                    ranges::begin(baseCoverages), std::plus<unsigned>{});
  cachedBaseWeights.clear();

  return *this;
}

RingmapData RingmapData::operator+(const RingmapData &other) const {
  RingmapData copy(*this);
  return copy += other;
}

void RingmapData::allowSequenceMismatches(bool value) {
#ifndef NDEBUG
  allowedMismatches = value;
#else
  (void)value;
#endif
}

void RingmapData::shuffle() { m_data.shuffle(); }

void RingmapData::resize(unsigned size) {
  // Fake coverage, should be fine for 'normal' operations
  ranges::transform(baseCoverages, ranges::begin(baseCoverages),
                    [oldSize = m_data.rows_size(), size](auto coverage) {
                      return static_cast<unsigned>(
                          static_cast<double>(coverage) / oldSize * size);
                    });

  if (size < nSetReads)
    nSetReads = size;
  m_data.resize(size);

  if (readsMap.size() != 0)
    readsMap.resize(size);
}

std::size_t RingmapData::size() const { return nSetReads; }

const std::vector<unsigned> &RingmapData::getBaseCoverages() const {
  return baseCoverages;
}

const std::vector<double> &RingmapData::getBaseWeights() const {
  assert(m_data.rows_size() > 0);
  assert(nSetReads == m_data.rows_size());

  if (cachedBaseWeights.size() > 0)
    return cachedBaseWeights;

  unsigned highestCoverage = *ranges::max_element(baseCoverages);
  cachedBaseWeights.resize(baseCoverages.size());
  ranges::transform(baseCoverages, ranges::begin(cachedBaseWeights),
                    [highestCoverage](unsigned coverage) {
                      return static_cast<double>(coverage) / highestCoverage;
                    });

  return cachedBaseWeights;
}

unsigned RingmapData::getModificationsFilter() const {
  return modificationsFilter;
}

bool RingmapData::operator==(const RingmapData &other) const {
  return startIndex == other.startIndex and endIndex == other.endIndex and
         sequence == other.sequence and basesFiltered == other.basesFiltered and
         modificationsFilter == other.modificationsFilter and
         nSetReads == other.nSetReads and
         minimumCoverage == other.minimumCoverage and
#ifndef NDEBUG
         allowedMismatches == other.allowedMismatches and
#endif
         m_data == other.m_data and oldColsToNew == other.oldColsToNew and
         readsMap == other.readsMap and baseCoverages == other.baseCoverages;
}

const std::string &RingmapData::getSequence() const { return sequence; }

WeightedClusters
RingmapData::getUnfilteredWeights(const WeightedClusters &weights) const {
  assert(weights.getElementsSize() == m_data.cols_size());

  if (basesFiltered) {
    WeightedClusters allWeights(sequence.size(), weights.getClustersSize(),
                                false);
    for (const auto &oldAndNewCol : oldColsToNew) {
      auto &&baseWeights = weights[oldAndNewCol.second];
      auto &&newBaseWeights = allWeights[oldAndNewCol.first];
      ranges::copy(baseWeights, ranges::begin(newBaseWeights));
    }

    return allWeights;
  } else
    return weights;
}

auto RingmapData::fractionReadsByWeights(const WeightedClusters &weights) const
    -> std::tuple<clusters_fraction_type, clusters_pattern_type,
                  clusters_assignment_type> {
  assert(weights.getElementsSize() == sequence.size());
  WeightedClusters reducedWeights;
  clusters_assignment_type clustersAssignment;
  const auto &usableWeights = [&, this]() -> decltype(auto) {
    if (basesFiltered) {
      reducedWeights = WeightedClusters(m_data.cols_size(),
                                        weights.getClustersSize(), false);
      for (const auto &oldAndNewCol : oldColsToNew) {
        assert(oldAndNewCol.first < weights.getElementsSize());
        assert(oldAndNewCol.second < reducedWeights.getElementsSize());

        auto &&weight = weights[oldAndNewCol.first];
        auto &&reducedWeight = reducedWeights[oldAndNewCol.second];
        assert(weight.span_size() == reducedWeight.span_size());

        ranges::copy(weight, ranges::begin(reducedWeight));
      }

      return const_cast<std::add_const_t<decltype(reducedWeights)> &>(
          reducedWeights);
    } else
      return (weights);
  }();

  struct rowHash {
    std::size_t operator()(
        std::reference_wrapper<const typename RingmapMatrix::row_type> rowRef)
        const {
      std::size_t hash = 0;
      std::size_t index = 0;
      for (auto &&value : rowRef.get())
        hash ^= std::hash<std::decay_t<decltype(value)>>{}(value) << index;

      return hash;
    }
  };

  struct rowEqual {
    using row_ref_type =
        std::reference_wrapper<const typename RingmapMatrix::row_type>;

    bool operator()(row_ref_type a, row_ref_type b) const {
      return a.get() == b.get();
    }
  };

  std::unordered_map<
      std::reference_wrapper<const typename RingmapMatrix::row_type>,
      std::size_t, rowHash, rowEqual>
      uniqueReadsCount;
  for (auto &&read : m_data.rows()) {
    auto indices = std::cref(read.modifiedIndices());
    if (auto uniqueReadIter = uniqueReadsCount.find(indices);
        uniqueReadIter != ranges::end(uniqueReadsCount))
      ++uniqueReadIter->second;
    else
      uniqueReadsCount[indices] = 1;
  }

  std::vector<std::size_t> counts(usableWeights.getClustersSize(), 0);
  std::vector<std::size_t> readsPerCluster(usableWeights.getClustersSize());
  clusters_pattern_type patterns(
      usableWeights.getClustersSize(),
      cluster_pattern_type(usableWeights.getElementsSize(), 0));
  {
    std::vector<double> readScore(usableWeights.getClustersSize());
    const auto readScoreEnd = ranges::end(readScore);
    std::vector<double> read(m_data.cols_size());
    for (auto &&pairedData : uniqueReadsCount) {
      auto &mutatedIndices = std::get<0>(pairedData);
      auto &&readOccurrences = std::get<1>(pairedData);
      ranges::fill(read, 0.);
      for (std::size_t index : mutatedIndices.get())
        read[index] = 1.;

      auto readScoreIter = ranges::begin(readScore);
      ranges::fill(readScoreIter, readScoreEnd, 0);
      for (auto &&cluster : usableWeights.clusters())
        *readScoreIter++ =
            ranges::inner_product(read, ranges::begin(cluster), 0.);

      auto assignedReadIter = [&]() {
        typename clusters_assignment_type::iterator assignmentIter =
            clustersAssignment.find(mutatedIndices);
        if (assignmentIter == ranges::end(clustersAssignment)) {
          assignmentIter =
              clustersAssignment
                  .emplace(mutatedIndices, std::vector<std::size_t>(
                                               weights.getClustersSize(), 0))
                  .first;
        }

        return assignmentIter;
      }();
      auto &assignedRead = assignedReadIter->second;

      std::vector<std::size_t> bestAssignmentIndices;
      {
        auto bestScoreIter = ranges::max_element(readScore);
        const auto bestScore = *bestScoreIter;
        for (; bestScoreIter != ranges::end(readScore);
             bestScoreIter = ranges::find(bestScoreIter, ranges::end(readScore),
                                          bestScore)) {
          bestAssignmentIndices.emplace_back(
              ranges::distance(ranges::begin(readScore), bestScoreIter++));
        }
      }

      const auto nAssignments =
          static_cast<unsigned>(bestAssignmentIndices.size());
      for (const std::size_t assignmentIndex : bestAssignmentIndices) {
        auto &pattern = patterns[assignmentIndex];

        ranges::transform(read, pattern, ranges::begin(pattern),
                          [readOccurrences = readOccurrences,
                           nAssignments](auto &&isMutated, std::size_t count) {
                            return static_cast<unsigned>(
                                static_cast<double>(count) +
                                isMutated *
                                    static_cast<double>(readOccurrences) /
                                    nAssignments);
                          });
        auto currentlyAssignedReads = readOccurrences / nAssignments;
        assignedRead[assignmentIndex] = currentlyAssignedReads;
        counts[assignmentIndex] += currentlyAssignedReads;
      }
    }
  }
  clusters_fraction_type fractions(weights.getClustersSize());
  auto const all_counts = ranges::accumulate(counts, std::size_t(0));
  if (all_counts != 0) {
    ranges::transform(
        counts, ranges::begin(fractions),
        [&, all_counts = static_cast<double>(all_counts)](std::size_t count) {
          return static_cast<double>(count) / all_counts;
        });
  } else if (fractions.size() == 1) {
    fractions[0] = 1.;
  }

  return std::tuple(std::move(fractions), std::move(patterns),
                    std::move(clustersAssignment));
}

void RingmapData::setBasesMask(arma::Col<std::uint8_t> mask) {
  basesMask = std::move(mask);
  if (not basesMask.empty() and basesMask.size() != sequence.size())
    throw std::runtime_error("mask must be of the same length of the sequence");
}

double RingmapData::getUnfoldedFraction() const {
  arma::mat frequencies = m_data.covariance();
  {
    arma::vec coverages(baseCoverages.size());
    ranges::copy(baseCoverages, ranges::begin(coverages));
    frequencies.each_col() /= coverages;
  }
  arma::mat unfoldedFrequencies = frequencies.diag() * frequencies.diag().t();
  assert(unfoldedFrequencies.size() == frequencies.size());

  struct InnerMatrix {
    InnerMatrix() = default;

    InnerMatrix(arma::mat &&plain) noexcept
        : plain(std::move(plain)), squared(arma::pow(this->plain, 2)) {}

    arma::mat plain;
    arma::mat squared;
  };

  unfoldedFrequencies.diag() = arma::zeros(frequencies.n_rows);
  frequencies.diag() = arma::zeros(frequencies.n_rows);

  double previousFraction = 0.1;
  double step = 0.01;

  auto calcInnerMatrix = [&](double fraction) {
    return InnerMatrix(frequencies - unfoldedFrequencies * fraction);
  };
  InnerMatrix previousInnerMatrix = calcInnerMatrix(previousFraction);

  auto functionAt = [](const InnerMatrix &innerMatrix) {
    return arma::accu(innerMatrix.squared);
  };

  auto calculateGradient = [&](const InnerMatrix &innerMatrix) {
    double gradient = arma::accu(innerMatrix.plain * unfoldedFrequencies * -2.);
    assert(not std::isnan(gradient));
    return gradient;
  };

  auto getDirection = [](double gradient) {
    return std::copysign(1., -gradient);
  };

  constexpr std::array<double, 2> wolfeParams{1e-4, 0.9};
  constexpr double epsilon = 1e-16;
  auto wolfeConditions = [&](double gradient, double step)
      -> std::optional<std::tuple<InnerMatrix, double, double>> {
    const double direction = getDirection(gradient);
    const double nextFraction = previousFraction - step * gradient;
    InnerMatrix nextInnerMatrix = calcInnerMatrix(nextFraction);
    if (functionAt(nextInnerMatrix) >
        functionAt(previousInnerMatrix) +
            wolfeParams[0] * step * direction * gradient)
      return {};

    const double nextGradient = calculateGradient(nextInnerMatrix);
    if (direction * nextGradient >= wolfeParams[1] * direction * gradient)
      return {std::make_tuple(std::move(nextInnerMatrix), nextFraction,
                              nextGradient)};
    else
      return {};
  };

  double previousGradient = calculateGradient(previousInnerMatrix);
  double currentGradient, currentFraction;
  InnerMatrix currentInnerMatrix;
  for (;;) {
    auto wolfeResult = wolfeConditions(previousGradient, step);
    if (wolfeResult) {
      std::tie(currentInnerMatrix, currentFraction, currentGradient) =
          std::move(*wolfeResult);
      break;
    }

    step *= 1.5;
    assert(not std::isinf(step));
  }

  auto calculateStep = [&] {
    return (currentFraction - previousFraction) /
           (currentGradient - previousGradient);
  };

  for (;;) {
    double currentValue = functionAt(currentInnerMatrix);
    double previousValue = functionAt(previousInnerMatrix);
    assert(currentValue < previousValue);

    if (previousValue - currentValue <= epsilon)
      break;

    step = calculateStep();
    previousFraction = std::exchange(currentFraction,
                                     currentFraction - step * currentGradient);
    previousInnerMatrix =
        std::exchange(currentInnerMatrix, calcInnerMatrix(currentFraction));
    previousGradient =
        std::exchange(currentGradient, calculateGradient(currentInnerMatrix));
  }

  return currentFraction;
}

auto RingmapData::remapPatterns(const clusters_pattern_type &patterns) const
    -> clusters_pattern_type {
  if (patterns.empty())
    return {};

  if (ranges::any_of(
          ranges::next(ranges::begin(patterns)), ranges::end(patterns),
          [&, frontSize = patterns.front().size()](const auto &pattern) {
            return pattern.size() != frontSize;
          }))
    throw std::logic_error("all patterns must have the same size");

  {
    std::size_t usableBases = 0;
    auto const sequence_size = static_cast<unsigned>(sequence.size());
    for (unsigned baseIndex = 0; baseIndex < sequence_size; ++baseIndex) {
      if ((basesMask.empty() or basesMask[baseIndex] == 1) and
          oldColsToNew.find(baseIndex) != ranges::end(oldColsToNew))
        ++usableBases;
    }

    if (patterns.front().size() != usableBases)
      throw std::logic_error(
          "patterns size is not related to this Ringmap data");
  }

  clusters_pattern_type mappedPatterns(
      patterns.size(),
      cluster_pattern_type(sequence.size(), static_cast<std::size_t>(0)));

  auto mappedPatternIter = ranges::begin(mappedPatterns);
  auto patternIter = ranges::begin(patterns);
  for (; patternIter != ranges::end(patterns);
       ++patternIter, ++mappedPatternIter) {
    const auto &pattern = *patternIter;
    auto &mappedPattern = *mappedPatternIter;

    for (unsigned baseIndex = 0,
                  sequence_size = static_cast<unsigned>(sequence.size());
         baseIndex < sequence_size; ++baseIndex) {
      if (not basesMask.empty() and basesMask[baseIndex] == 0.)
        continue;

      if (auto newColIter = oldColsToNew.find(baseIndex);
          newColIter != ranges::end(oldColsToNew))
        mappedPattern[baseIndex] = pattern[newColIter->second];
    }
  }

  return mappedPatterns;
}

void RingmapData::enqueueRingmapsFromMutationMap(
    MutationMap &mutationMap,
    parallel::blocking_queue<std::pair<MutationMapTranscript, RingmapData>>
        &queue,
    Args const &args) {

  if (auto const &whitelist_filename = args.whitelist();
      whitelist_filename.empty()) {
    for (auto &&transcript : mutationMap)
      queue.push(std::pair(transcript, RingmapData(transcript, args)));
  } else {
    std::vector<std::string> whitelisted_genes;
    {
      std::ifstream whitelist_stream(whitelist_filename);
      for (std::string line; std::getline(whitelist_stream, line);
           line.clear()) {
        ranges::transform(line, ranges::begin(line), [](char c) {
          return static_cast<char>(std::tolower(c));
        });
        whitelisted_genes.emplace_back(std::move(line));
      }
    }

    ranges::sort(whitelisted_genes);
    std::string transcriptId;
    for (auto &&transcript : mutationMap) {
      transcriptId = transcript.getId();
      ranges::transform(transcriptId, ranges::begin(transcriptId), [](char c) {
        return static_cast<char>(std::tolower(c));
      });

      if (ranges::binary_search(whitelisted_genes, transcriptId)) {
        queue.push(std::pair(transcript, RingmapData(transcript, args)));
      }
    }
  }

  queue.finish();
}

RingmapData RingmapData::get_new_range(
    unsigned begin, unsigned end,
    std::vector<unsigned> *const used_reads_indices) const {
  assert(not basesFiltered);
  assert(oldColsToNew.empty());
  assert(begin >= startIndex);
  assert(end <= endIndex);
  assert(basesMask.empty());
  assert(readsMap.empty());

  RingmapData new_ringmap;
#ifndef NDEBUG
  new_ringmap.allowedMismatches = allowedMismatches;
#endif
  new_ringmap.startIndex = begin;
  new_ringmap.endIndex = end;
  new_ringmap.modificationsFilter = modificationsFilter;
  new_ringmap.minimumCoverage = minimumCoverage;
  new_ringmap.minimumModificationsPerBase = minimumModificationsPerBase;
  new_ringmap.minimumModificationsPerRead = minimumModificationsPerRead;
  new_ringmap.minimumModificationsPerBaseFraction =
      minimumModificationsPerBaseFraction;
  new_ringmap.sequence = sequence.substr(begin - startIndex, end - begin);
  new_ringmap.baseCoverages.resize(end - begin, 0);
  new_ringmap.m_data = RingmapMatrix(end - begin);
  new_ringmap.shape = shape;

  auto &&rows = m_data.rows();
  auto rows_iter = ranges::begin(rows);
  auto const end_rows = ranges::end(rows);

  for (unsigned row_index = 0; rows_iter < end_rows; ++rows_iter, ++row_index) {
    auto &&row = *rows_iter;

    assert(row.begin_index() <= row.end_index());
    if (row.begin_index() <= begin and row.end_index() >= end) {
      if (used_reads_indices) {
        used_reads_indices->emplace_back(row_index);
      }

      auto &&indices = row.modifiedIndices();
      assert(ranges::is_sorted(indices));

      auto const first_index_iter = ranges::lower_bound(indices, begin);
      auto const last_index_iter =
          ranges::lower_bound(first_index_iter, ranges::end(indices), end);

      ringmap_matrix::row_type new_row(static_cast<std::size_t>(
          ranges::distance(first_index_iter, last_index_iter)));
      new_row.begin_index = begin;
      new_row.end_index = end;

      assert(
          ranges::all_of(first_index_iter, last_index_iter, [&](auto &&index) {
            return index >= begin and index < end;
          }));
      ranges::transform(first_index_iter, last_index_iter,
                        ranges::begin(new_row),
                        [begin](auto &&index) { return index - begin; });
      new_ringmap.m_data.addModifiedIndicesRow(std::move(new_row));
    }
  }

  new_ringmap.nSetReads = static_cast<unsigned>(new_ringmap.m_data.rows_size());
  ranges::fill(new_ringmap.baseCoverages, new_ringmap.nSetReads);

  return new_ringmap;
}

std::vector<RingmapData> RingmapData::split_into_windows(
    std::vector<results::Window> const &windows) && {
  assert(not windows.empty());
  auto const windows_size = windows.size();
  std::vector<RingmapData> splitted_ringmaps(windows_size);

  if (windows_size == 1 and windows[0].begin_index == startIndex and
      windows[0].end_index == endIndex) {

    splitted_ringmaps[0] = std::move(*this);
    return splitted_ringmaps;
  }

  auto const windows_end = std::end(windows);
  {
    auto windows_iter = std::begin(windows);
    auto splitted_ringmaps_iter = std::begin(splitted_ringmaps);
    for (; windows_iter < windows_end;
         ++windows_iter, ++splitted_ringmaps_iter) {
      auto &window = *windows_iter;
      auto &ringmap = *splitted_ringmaps_iter;

      auto const window_size =
          static_cast<unsigned>(window.end_index - window.begin_index);
#ifndef NDEBUG
      ringmap.allowedMismatches = allowedMismatches;
#endif
      ringmap.startIndex = window.begin_index;
      ringmap.endIndex = window.end_index;
      ringmap.modificationsFilter = modificationsFilter;
      ringmap.minimumCoverage = minimumCoverage;
      ringmap.minimumModificationsPerBase = minimumModificationsPerBase;
      ringmap.sequence =
          sequence.substr(window.begin_index - startIndex, window_size);
      ringmap.baseCoverages.resize(window_size, 0);
      ringmap.m_data = RingmapMatrix(window_size);
    }
  }

  auto &&rows = m_data.rows();
  auto rows_iter = ranges::begin(rows);
  auto const end_rows = ranges::end(rows);

  auto const calc_modifications_intersection = [](auto const &row,
                                                  auto const &window) {
    if (row.end_index() <= window.begin_index or
        row.begin_index() >= window.end_index)
      return 0u;

    return static_cast<unsigned>(
        ranges::count_if(row.modifiedIndices(), [&](auto base_index) {
          return base_index >= window.begin_index and
                 base_index < window.end_index;
        }));
  };

  auto const calc_intersection = [](auto const &row,
                                    auto const &window) -> unsigned {
    auto const row_begin = row.begin_index();
    auto const row_end = row.end_index();

    if (row_begin <= window.begin_index) {
      if (row_end <= window.end_index) {
        return row_end - window.begin_index;
      } else {
        return window.end_index - window.begin_index;
      }
    } else if (row_begin < window.end_index) {
      if (row_end <= window.end_index) {
        return row_end - row_begin;
      } else {
        return window.end_index - row_begin;
      }
    } else {
      return 0u;
    }
  };

  std::vector<std::size_t> candidate_indices, last_step_indices;
  for (; rows_iter < end_rows; ++rows_iter) {
    candidate_indices.clear();
    auto &&row = *rows_iter;

    {
      // Put the first window as a candidate
      candidate_indices.resize(1);
      candidate_indices[0] = 0;
      unsigned best_intersection_size =
          calc_modifications_intersection(row, windows[0]);

      for (std::size_t window_index = 1; window_index < windows_size;
           ++window_index) {
        if (auto intersection_value =
                calc_modifications_intersection(row, windows[window_index]);
            intersection_value > best_intersection_size) {
          candidate_indices.resize(1);
          candidate_indices[0] = window_index;
          best_intersection_size = intersection_value;
        } else if (intersection_value == best_intersection_size) {
          candidate_indices.emplace_back(window_index);
        }
      }

      if (best_intersection_size == 0u)
        continue;
    }

    // More than one candidate using modified bases. Using read intersection.
    if (candidate_indices.size() > 1) {
      std::swap(candidate_indices, last_step_indices);
      candidate_indices.resize(1);
      candidate_indices[0] = last_step_indices[0];
      unsigned best_intersection_size =
          calc_intersection(row, windows[last_step_indices[0]]);

      std::for_each(std::next(std::begin(last_step_indices)),
                    std::end(last_step_indices), [&](auto window_index) {
                      if (auto intersection_value =
                              calc_intersection(row, windows[window_index]);
                          intersection_value > best_intersection_size) {
                        candidate_indices.resize(1);
                        candidate_indices[0] = window_index;
                      } else if (intersection_value == best_intersection_size) {
                        candidate_indices.emplace_back(window_index);
                      }
                    });
    }

    // Still more than one candidate? Let's take the smaller window
    if (candidate_indices.size() > 1) {
      auto const best_candidate_index =
          *ranges::min_element(candidate_indices, {}, [&](auto index) {
            auto const &window = windows[index];
            return window.end_index - window.begin_index;
          });

      candidate_indices[0] = best_candidate_index;
    }

    auto const best_candidate_index = candidate_indices[0];
    auto &ringmap = splitted_ringmaps[best_candidate_index];
    auto const &window = windows[best_candidate_index];

    auto const &indices = row.modifiedIndices();
    assert(ranges::is_sorted(indices));

    auto const new_row_begin =
        std::max(row.begin_index(), static_cast<unsigned>(window.begin_index));
    auto const new_row_end =
        std::min(row.end_index(), static_cast<unsigned>(window.end_index));

    auto const first_index_iter = ranges::lower_bound(indices, new_row_begin);
    auto const last_index_iter = ranges::lower_bound(
        first_index_iter, ranges::end(indices), new_row_end);

    ringmap_matrix::row_type new_row(static_cast<std::size_t>(
        ranges::distance(first_index_iter, last_index_iter)));
    new_row.begin_index = new_row_begin;
    new_row.end_index = new_row_end;

    assert(ranges::all_of(first_index_iter, last_index_iter, [&](auto &&index) {
      return index >= new_row_begin and index < new_row_end;
    }));
    ranges::transform(
        first_index_iter, last_index_iter, ranges::begin(new_row),
        [begin = window.begin_index](auto &&index) { return index - begin; });
    ringmap.m_data.addModifiedIndicesRow(std::move(new_row));

    auto const base_coverages_begin = std::begin(ringmap.baseCoverages);
    auto const first_base_cov =
        std::next(base_coverages_begin, new_row_begin - window.begin_index);
    std::transform(
        first_base_cov,
        std::next(base_coverages_begin, new_row_end - window.begin_index),
        first_base_cov, [](auto base_coverage) { return base_coverage + 1; });
  }

  for (auto &ringmap : splitted_ringmaps) {
    ringmap.nSetReads = static_cast<unsigned>(ringmap.m_data.rows_size());
  }

  return splitted_ringmaps;
}

#pragma once

#include "args.hpp"
#include "parallel/blocking_queue.hpp"
#include "ringmap_matrix.hpp"
#include "weighted_clusters_impl.hpp"

#include <array>
#include <map>
#include <string>
#include <vector>

class MutationMap;
class MutationMapTranscript;
class RnaSecondaryStructure;
struct PairedRnaSecondaryStructure;

namespace test {
struct RingmapData;
} // namespace test

namespace results {
struct Window;
} // namespace results

class RingmapData {
public:
  using data_type = RingmapMatrix;
  using data_value_type = RingmapMatrix::value_type;
  friend struct test::RingmapData;

  RingmapData() = default;
  RingmapData(const std::string& sequence, data_type&& dataMatrix,
              unsigned startIndex, unsigned endIndex, Args const& args);
  RingmapData(const std::string& filename, const std::string& sequence,
              Args const& args, bool keepFragments = true);
  RingmapData(const MutationMapTranscript& transcript, Args const& args);
  template <typename Iter>
  RingmapData(const std::string& sequence, unsigned nReads, Iter readsBegin,
              Iter readsEnd, Args const& args);
  RingmapData(const RingmapData& other,
              const std::vector<unsigned>& subsetIndices);

  template <typename T>
  static std::enable_if_t<std::is_same<std::decay_t<T>, RingmapData>::value,
                          RingmapData>
  from(T&& other);

  bool isFiltered() const;
  unsigned getModificationsFilter() const;
  void filterBases();
  void filterReads();
  void filter();
  void perturb();
  void shuffle();
  void resize(unsigned size);
  const data_type& data() const;
  static void removeHighValuesOnAdjacency(arma::mat& adjacency,
                                          double maxFraction = 0.35);
  void shed(unsigned begin,
            unsigned end = std::numeric_limits<unsigned>::max());
  void remove(unsigned begin,
              unsigned end = std::numeric_limits<unsigned>::max());
  void dump(const std::string& filename) const;
  std::size_t size() const;
  const std::string& getSequence() const;

  void allowSequenceMismatches(bool value = true);
  RingmapData& operator+=(const RingmapData& other);
  RingmapData operator+(const RingmapData& other) const;
  bool operator==(const RingmapData& other) const;

  template <typename Iterable>
  void keepOnlyIndices(Iterable&& iterable);

  void setBasesMask(arma::Col<std::uint8_t> mask);
  RingmapData get_new_range(unsigned begin, unsigned end,
                            std::vector<unsigned>* const = nullptr) const;

  WeightedClusters getUnfilteredWeights(const WeightedClusters& weights) const;

  using clusters_fraction_type = std::vector<double>;
  using cluster_pattern_type = std::vector<std::size_t>;
  using clusters_pattern_type = std::vector<cluster_pattern_type>;
  using clusters_assignment_type =
      std::map<data_type::row_type, std::vector<std::size_t>>;
  std::tuple<clusters_fraction_type, clusters_pattern_type,
             clusters_assignment_type>
  fractionReadsByWeights(const WeightedClusters& weights) const;
  clusters_pattern_type
  remapPatterns(const clusters_pattern_type& patterns) const;

  double getUnfoldedFraction() const;

  static void enqueueRingmapsFromMutationMap(
      MutationMap& mutationMap,
      parallel::blocking_queue<std::pair<MutationMapTranscript, RingmapData>>&
          queue,
      Args const& args);

private:
  RingmapData(const std::string& sequence, data_type&& dataMatrix,
              unsigned startIndex, unsigned endIndex);

#ifndef NDEBUG
  bool allowedMismatches = false;
#endif

  bool basesFiltered = false;
  unsigned startIndex, endIndex;
  unsigned nSetReads = 0;
  unsigned modificationsFilter = 0;
  unsigned minimumCoverage = 100;
  unsigned minimumModificationsPerBase = 2;
  unsigned minimumModificationsPerRead = 2;
  float minimumModificationsPerBaseFraction = 0.005f;
  std::string sequence;
  std::vector<unsigned> readsMap;
  std::vector<unsigned> baseCoverages;
  mutable std::vector<double> cachedBaseWeights;
  std::map<unsigned, unsigned> oldColsToNew;
  data_type m_data;
  arma::Col<std::uint8_t> basesMask;
  bool shape = false;

public:
  const decltype(oldColsToNew)& getNonFilteredToFilteredMap() const;
  decltype(oldColsToNew) getFilteredToNonFilteredMap() const;
  const decltype(readsMap)& getReadsMap() const;
  const std::vector<unsigned>& getBaseCoverages() const;
  const std::vector<double>& getBaseWeights() const;
  std::vector<RingmapData>
  split_into_windows(std::vector<results::Window> const& windows) &&;

private:
  template <typename Iter>
  void addReads(Iter begin, Iter end);
};

#include "ringmap_data_impl.hpp"
#include "test/ringmap_data.hpp"

#pragma once

#include "args.hpp"
#include "concepts.hpp"
#include "parallel/blocking_queue.hpp"
#include "read_clusters_assignments.hpp"
#include "ringmap_matrix.hpp"
#include "weighted_clusters_impl.hpp"

#include <boost/thread/lock_options.hpp>
#include <boost/thread/lock_types.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
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

struct CachedBaseWeights {
  using data_type = std::vector<double>;

  inline CachedBaseWeights() noexcept = default;
  inline CachedBaseWeights(CachedBaseWeights const &other)
      : weights_(([&] {
          boost::shared_lock lock(other.mutex_);
          if (other.weights_) {
            return std::make_shared<data_type>(*other.weights_);
          } else {
            return std::shared_ptr<data_type>();
          }
        })()) {}
  inline CachedBaseWeights(CachedBaseWeights &&other) noexcept
      : weights_(([&] {
          static_assert(std::is_nothrow_move_constructible_v<
                        std::shared_ptr<std::vector<double>>>);

          boost::unique_lock lock(other.mutex_);
          if (other.weights_) {
            return std::shared_ptr(std::move(other.weights_));
          } else {
            return std::shared_ptr<data_type>();
          }
        })()) {}

  inline CachedBaseWeights &operator=(CachedBaseWeights const &other) {
    if (this == &other) {
      return *this;
    }

    boost::unique_lock this_lock(mutex_, boost::defer_lock);
    boost::shared_lock other_lock(other.mutex_, boost::defer_lock);
    if (&this->mutex_ < &other.mutex_) {
      this_lock.lock();
      other_lock.lock();
    } else {
      other_lock.lock();
      this_lock.lock();
    }
    if (other.weights_) {
      weights_ = std::make_shared<std::vector<double>>(*other.weights_);
    } else {
      weights_.reset();
    }
    return *this;
  }

  inline CachedBaseWeights &operator=(CachedBaseWeights &&other) {
    static_assert(std::is_nothrow_move_assignable_v<
                  std::shared_ptr<std::vector<double>>>);
    if (this == &other) {
      return *this;
    }

    std::scoped_lock lock(mutex_, other.mutex_);
    if (other.weights_) {
      weights_ = std::move(other.weights_);
    } else {
      weights_.reset();
    }
    return *this;
  }

  constexpr void reset() const noexcept {
    boost::unique_lock lock(mutex_);
    weights_.reset();
  }

  template <typename F>
    requires InvocableR<F, data_type>
  constexpr std::weak_ptr<data_type> get_weak_or_init(F &&f) const {
    boost::upgrade_lock upgrade_lock(mutex_);
    if (weights_) {
      return std::weak_ptr(weights_);
    } else {
      boost::unique_lock unique_lock(std::move(upgrade_lock));
      weights_ = std::make_shared<data_type>(f());
      return std::weak_ptr(weights_);
    }
  }

protected:
  mutable boost::upgrade_mutex mutex_;
  mutable std::shared_ptr<data_type> weights_;
};

class RingmapData {
public:
  using data_type = RingmapMatrix;
  using data_value_type = RingmapMatrix::value_type;
  friend struct test::RingmapData;

  RingmapData() = default;
  RingmapData(std::string_view sequence, data_type &&dataMatrix,
              unsigned startIndex, unsigned endIndex, Args const &args);
  RingmapData(const std::string &filename, std::string_view sequence,
              Args const &args, bool keepFragments = true);
  RingmapData(const MutationMapTranscript &transcript, Args const &args);
  template <typename Iter>
  RingmapData(std::string_view sequence, unsigned nReads, Iter readsBegin,
              Iter readsEnd, Args const &args);
  RingmapData(const RingmapData &other,
              const std::vector<unsigned> &subsetIndices);

  template <typename T>
  static std::enable_if_t<std::is_same<std::decay_t<T>, RingmapData>::value,
                          RingmapData>
  from(T &&other);

  bool isFiltered() const;
  unsigned getModificationsFilter() const;
  void filterBases();
  void filterReads();
  void filter();
  void perturb();
  void shuffle();
  void resize(unsigned size);
  const data_type &data() const;
  static void removeHighValuesOnAdjacency(arma::mat &adjacency,
                                          double maxFraction);
  void shed(unsigned begin,
            unsigned end = std::numeric_limits<unsigned>::max());
  void remove(unsigned begin,
              unsigned end = std::numeric_limits<unsigned>::max());
  void dump(const std::string &filename) const;
  std::size_t size() const;
  const std::string &getSequence() const;

  void allowSequenceMismatches(bool value = true);
  RingmapData &operator+=(const RingmapData &other);
  RingmapData operator+(const RingmapData &other) const;
  bool operator==(const RingmapData &other) const;

  template <typename Iterable> void keepOnlyIndices(Iterable &&iterable);

  void setBasesMask(arma::Col<std::uint8_t> mask);
  RingmapData get_new_range(unsigned begin, unsigned end,
                            std::vector<unsigned> *const = nullptr) const;

  WeightedClusters getUnfilteredWeights(const WeightedClusters &weights) const;

  using clusters_fraction_type = std::vector<double>;
  using cluster_pattern_type = std::vector<std::size_t>;
  using clusters_pattern_type = std::vector<cluster_pattern_type>;
  using clusters_assignment_type =
      std::map<data_type::row_type, ReadClustersAssignments>;
  std::tuple<clusters_fraction_type, clusters_pattern_type,
             clusters_assignment_type>
  fractionReadsByWeights(const WeightedClusters &weights, unsigned window_size,
                         bool skip_ambiguous_assignments) const;
  clusters_pattern_type
  remapPatterns(const clusters_pattern_type &patterns) const;

  double getUnfoldedFraction() const;

  static void enqueueRingmapsFromMutationMap(
      MutationMap &mutationMap,
      parallel::blocking_queue<std::tuple<MutationMapTranscript, RingmapData,
                                          std::optional<std::uint16_t>>> &queue,
      Args const &args);

  template <typename R>
    requires std::ranges::range<R> and
             std::same_as<std::ranges::range_reference_t<R>, RingmapData &> and
             std::random_access_iterator<std::ranges::iterator_t<R>>
  static void filter_bases_on_replicates(R &ringmap_data_range);

  template <typename R>
    requires std::ranges::range<R> and
             std::same_as<std::ranges::range_reference_t<R>, RingmapData &> and
             std::random_access_iterator<std::ranges::iterator_t<R>>
  static void filter_bases_on_replicates_for_assignments(R &ringmap_data_range);

  constexpr Args const &args() const noexcept { return *args_; }

private:
  RingmapData(std::string_view sequence, data_type &&dataMatrix,
              unsigned startIndex, unsigned endIndex);

  template <typename R, typename F>
    requires std::ranges::range<R> and
             std::same_as<std::ranges::range_reference_t<R>, RingmapData &> and
             std::random_access_iterator<std::ranges::iterator_t<R>>
  static void filter_bases_on_replicates_inner(R &ringmap_data_range,
                                               F base_filter);

#ifndef NDEBUG
  bool allowedMismatches = false;
#endif

  bool basesFiltered = false;
  unsigned startIndex, endIndex;
  unsigned nSetReads = 0;
  unsigned modificationsFilter = 0;
  std::string sequence;
  std::vector<unsigned> readsMap;
  std::vector<unsigned> baseCoverages;
  CachedBaseWeights cached_base_weights;

  std::map<unsigned, unsigned> oldColsToNew;
  data_type m_data;
  arma::Col<std::uint8_t> basesMask;
  Args const *args_;

public:
  const decltype(oldColsToNew) &getNonFilteredToFilteredMap() const;
  decltype(oldColsToNew) getFilteredToNonFilteredMap() const;
  const decltype(readsMap) &getReadsMap() const;
  const std::vector<unsigned> &getBaseCoverages() const;
  const std::weak_ptr<std::vector<double>> getBaseWeights() const;
  std::vector<RingmapData>
  split_into_windows(std::vector<results::Window> const &windows) &&;

private:
  template <typename Iter> void addReads(Iter begin, Iter end);
};

#include "ringmap_data_impl.hpp"
#include "test/ringmap_data.hpp"

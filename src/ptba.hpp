#pragma once

#include "args.hpp"
#include "draco.hpp"
#include "ptba_types.hpp"
#include "results/transcript.hpp"
#include "ringmap_data.hpp"

#include <cmath>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include <armadillo>

namespace log_data {

struct NotEnoughReads {
  std::size_t reads;
  std::size_t min_filtered_reads;
};

struct NotEnoughBases {
  std::size_t bases;
  std::size_t min_bases;
};

struct AllZeroEigenGaps {};

struct NoUsefulEigenGaps {
  std::size_t eigengaps;
  std::size_t max_clusters;
};

namespace permuting {

struct NotEnoughReads {
  unsigned permutation;
  std::size_t filtered_reads;
  std::size_t min_filtered_reads;
};

struct NewValidEigengap {
  unsigned permutation;
  std::optional<std::size_t> previous_eigengap_index;
  std::size_t new_eigengap_index;
  bool distribution_is_evaluated;
  double eigengap_difference;
  double cumulative_difference;
  double min_eigengap_threshold;
};

struct SignificantEigengapLowDifference {
  unsigned permutation;
  double eigengap_difference;
  double eigengap_diff_absolute_threshold;
  bool distribution_is_evaluated;
  std::optional<std::size_t> valid_eigengap_index;
  std::size_t current_eigengap_index;
};

struct NotSignificantEigengap {
  unsigned permutation;
  std::optional<std::size_t> valid_eigengap_index;
  std::size_t current_eigengap_index;
};

namespace solution_found {

struct UsefulEigengapIndex {};

struct OverExtendedSearch {
  std::size_t valid_eigengap_index;
  std::size_t extended_search_eigengaps;
};

using Solution = std::variant<solution_found::UsefulEigengapIndex,
                              solution_found::OverExtendedSearch>;

} // namespace solution_found

struct SolutionFound {
  unsigned permutation;
  std::size_t eigengap_index;
  std::size_t useful_eigengaps;
  solution_found::Solution solution;
};

using Data = std::variant<NotEnoughReads, NewValidEigengap,
                          SignificantEigengapLowDifference,
                          NotSignificantEigengap, SolutionFound>;
} // namespace permuting

struct Permuting {
  std::vector<permuting::Data> data;
};

} // namespace log_data

using LogData = std::variant<log_data::NotEnoughReads, log_data::NotEnoughBases,
                             log_data::AllZeroEigenGaps,
                             log_data::NoUsefulEigenGaps, log_data::Permuting>;

struct PtbaResult {
  unsigned nClusters;
  arma::mat eigenVecs;
  arma::vec eigenGaps;
  PerturbedEigengaps perturbedEigenGaps;
  std::vector<unsigned> significantIndices;
  std::vector<unsigned> filteredToUnfilteredBases;
  arma::mat adjacency;
  LogData log_data;

  PtbaResult(unsigned nClusters, arma::mat eigenVecs, arma::vec eigenGaps,
             PerturbedEigengaps perturbedEigenGaps,
             std::vector<unsigned> significantIndices,
             std::vector<unsigned> filteredToUnfilteredBases,
             arma::mat adjacency, LogData log_data)
      : nClusters(nClusters), eigenVecs(std::move(eigenVecs)),
        eigenGaps(std::move(eigenGaps)),
        perturbedEigenGaps(std::move(perturbedEigenGaps)),
        significantIndices(std::move(significantIndices)),
        filteredToUnfilteredBases(std::move(filteredToUnfilteredBases)),
        adjacency(std::move(adjacency)), log_data(std::move(log_data)) {}

  PtbaResult(LogData log_data) : log_data(std::move(log_data)) {}
  PtbaResult(const PtbaResult &) = default;
  PtbaResult(PtbaResult &&) = default;
  PtbaResult &operator=(const PtbaResult &) = default;
  PtbaResult &operator=(PtbaResult &&) = default;
};

class Ptba /* Permutation test-based analysis */
{
public:
  struct exception : std::runtime_error {
    explicit exception(const std::string &what) : std::runtime_error(what) {};
    explicit exception(const char *what) : std::runtime_error(what) {};
  };

  Ptba(const RingmapData &data, Args const &args);

  unsigned getNumberOfClustersAndDumpData(
      std::string_view eigenVecsFilename, std::string_view eigenGapsFilename,
      std::string_view perturbedEigenGapsFilename) const;
  std::pair<unsigned, std::vector<unsigned>>
  getNumberOfClustersAndSignificantIndices() const;
  std::pair<unsigned, std::vector<unsigned>>
  getNumberOfClustersAndSignificantIndicesAndDumpData(
      std::string_view eigenValsFilename,
      std::string_view perturbedEigenValsFilename) const;
  std::vector<unsigned> getAllSignificantEigenGapIndices() const;

  PtbaResult run() const noexcept(false);

  static void dumpEigenVecs(const arma::mat &eigenVecs,
                            std::string_view eigenVecsFilename);
  static void dumpEigenGaps(const arma::vec &eigenGaps,
                            std::string_view eigenGapsFilename);
  static void
  dumpPerturbedEigenGaps(PerturbedEigengaps const &perturbed_eigengaps,
                         std::string_view perturbedEigenGapsFilename);

private:
  static std::tuple<arma::mat, arma::vec, arma::vec, arma::mat>
  calculateEigenGaps(const RingmapData &data);

  template <typename Distribution>
  static bool
  is_distribution(Distribution &&distribution,
                  PerturbedEigengap const &perturbed_data) noexcept(false);

  const RingmapData *ringmapData;
  unsigned minFilteredReads = 5;
  unsigned maxPermutations = 400;
  unsigned minPermutations = 8;
  double minEigenGapThreshold = 0.10;
  double eigenGapDiffAbsoluteThreshold = 0.03;
  double alphaValue = 0.01;
  double betaValue = 0.2;
  unsigned maxClusters = std::numeric_limits<unsigned>::max();
  unsigned alternative_check_permutations = 50;
  double min_null_stddev = 0.025;
  unsigned minBasesSize = 10;
  unsigned char extended_search_eigengaps;
  bool ignore_first_eigengap;
};

void print_log_data(LogData const &log_data, Window const &window,
                    size_t window_index, unsigned window_size,
                    results::Transcript const &transcript,
                    std::size_t replicate_index) noexcept;

#include "ptba_impl.hpp"

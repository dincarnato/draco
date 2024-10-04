#pragma once

#include "args.hpp"
#include "ptba_types.hpp"
#include "ringmap_data.hpp"

#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <armadillo>

struct PtbaResult {
  unsigned nClusters;
  arma::mat eigenVecs;
  arma::vec eigenGaps;
  PerturbedEigengaps perturbedEigenGaps;
  std::vector<unsigned> significantIndices;
  std::vector<unsigned> filteredToUnfilteredBases;
  arma::mat adjacency;
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
  void setMaxClusters(unsigned value);
  std::vector<unsigned> getAllSignificantEigenGapIndices() const;
  void setMinEigenGapThreshold(double value);

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
  double firstEigengapThreshold = 0.95;
  double minEigenGapThreshold = 0.10;
  double eigenGapDiffAbsoluteThreshold = 0.03;
  double alphaValue = 0.01;
  double betaValue = 0.2;
  double firstEigengapBetaValue = 0.4;
  unsigned maxClusters = std::numeric_limits<unsigned>::max();
  unsigned alternative_check_permutations = 50;
  double min_null_stddev = 0.025;
  unsigned minBasesSize = 10;
  unsigned char extended_search_eigengaps;
};

#include "ptba_impl.hpp"

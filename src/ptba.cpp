#include "ptba.hpp"
#include "weibull_fitter.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iterator>

#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#error "Missing filesystem header"
#endif

Ptba::Ptba(const RingmapData &data, Args const &args)
    : ringmapData(&data), minFilteredReads(args.min_filtered_reads()),
      maxPermutations(args.max_permutations()),
      minPermutations(args.min_permutations()),
      firstEigengapThreshold(args.first_eigengap_threshold()),
      minEigenGapThreshold(args.min_eigengap_threshold()),
      eigenGapDiffAbsoluteThreshold(args.eigengap_diff_absolute_threshold()),
      alphaValue(args.alpha_value()), betaValue(args.beta_value()),
      firstEigengapBetaValue(args.first_eigengap_beta_value()),
      maxClusters(args.max_clusters()),
      alternative_check_permutations(args.alternative_check_permutations()),
      min_null_stddev(args.min_null_stddev()),
      minBasesSize(args.min_bases_size()),
      extended_search_eigengaps(args.extended_search_eigengaps()) {}

std::tuple<arma::mat, arma::vec, arma::vec, arma::mat>
Ptba::calculateEigenGaps(const RingmapData &data) {
  arma::mat normalizedLaplacian;
  arma::mat adjacency = data.data().covariance(data.getBaseWeights());
  {
    // data.fixBadNeighboursOnAdjacency(adjacency);
    RingmapData::removeHighValuesOnAdjacency(adjacency);

    adjacency.diag() = arma::zeros(adjacency.n_cols);
    arma::vec degree = arma::sum(adjacency, 1);
    arma::mat laplacian = arma::diagmat(degree) - adjacency;
    std::transform(std::begin(degree), std::end(degree), std::begin(degree),
                   [](double degree) {
                     if (std::abs(degree) < 1e-4)
                       return 1.;
                     else
                       return degree;
                   });
    arma::mat invSqrtDiagonal =
        arma::diagmat(arma::ones(degree.size()) / arma::sqrt(degree));
    normalizedLaplacian = invSqrtDiagonal * laplacian * invSqrtDiagonal;

    assert(not normalizedLaplacian.has_nan());
  }

  arma::mat eigVecs;
  arma::vec eigValues;
  arma::eig_sym(eigValues, eigVecs, normalizedLaplacian);
  return std::make_tuple(std::move(eigVecs), eigValues, arma::diff(eigValues),
                         std::move(adjacency));
}

unsigned Ptba::getNumberOfClusters() const { return run(); }

void Ptba::dumpEigenVecs(const arma::mat &eigenVecs,
                         std::string_view eigenVecsFilename) {
  std::ofstream eigenVecsStream{fs::path{eigenVecsFilename}};
  eigenVecsStream << std::setprecision(10) << std::scientific;
  for (std::size_t col = 0; col < eigenVecs.n_cols; ++col) {
    for (std::size_t row = 0; row < eigenVecs.n_rows - 1; ++row)
      eigenVecsStream << eigenVecs(row, col) << ' ';

    eigenVecsStream << eigenVecs(eigenVecs.n_rows - 1, col) << '\n';
  }
}

void Ptba::dumpEigenGaps(const arma::vec &eigenGaps,
                         std::string_view eigenGapsFilename) {
  std::ofstream eigenGapsStream{fs::path{eigenGapsFilename}};
  eigenGapsStream << std::setprecision(10) << std::scientific;
  unsigned dataIndex = 1;
  for (double eigenGap : eigenGaps)
    eigenGapsStream << dataIndex++ << ' ' << eigenGap << '\n';
}

void Ptba::dumpPerturbedEigenGaps(PerturbedEigengaps const &perturbed_eigengaps,
                                  std::string_view perturbedEigenGapsFilename) {
  std::ofstream eigenGapsStream{fs::path{perturbedEigenGapsFilename}};
  eigenGapsStream << std::setprecision(10) << std::scientific;
  unsigned eigengap_index = 1;
  for (auto &&perturbed_eigengap : perturbed_eigengaps) {
    for (auto &&value : perturbed_eigengap)
      eigenGapsStream << eigengap_index << ' ' << value << '\n';

    ++eigengap_index;
  }
}

unsigned Ptba::getNumberOfClustersAndDumpData(
    std::string_view eigenVecsFilename, std::string_view eigenGapsFilename,
    std::string_view perturbedEigenGapsFilename) const {
  PtbaResult result = result_from_run();
  dumpEigenVecs(result.eigenVecs, std::move(eigenVecsFilename));
  dumpEigenGaps(result.eigenGaps, std::move(eigenGapsFilename));
  dumpPerturbedEigenGaps(result.perturbedEigenGaps,
                         std::move(perturbedEigenGapsFilename));

  return result.nClusters;
}

std::pair<unsigned, std::vector<unsigned>>
Ptba::getNumberOfClustersAndSignificantIndices() const {
  PtbaResult result = result_from_run();
  return {result.nClusters, std::move(result.significantIndices)};
}

std::pair<unsigned, std::vector<unsigned>>
Ptba::getNumberOfClustersAndSignificantIndicesAndDumpData(
    std::string_view eigenGapsFilename,
    std::string_view perturbedEigenGapsFilename) const {
  PtbaResult result = result_from_run();
  dumpEigenGaps(result.eigenGaps, std::move(eigenGapsFilename));
  dumpPerturbedEigenGaps(result.perturbedEigenGaps,
                         std::move(perturbedEigenGapsFilename));

  return {result.nClusters, std::move(result.significantIndices)};
}

void Ptba::setMaxClusters(unsigned value) { maxClusters = value; }

std::vector<unsigned> Ptba::getAllSignificantEigenGapIndices() const {
  return result_from_run().significantIndices;
}

void Ptba::setMinEigenGapThreshold(double value) {
  minEigenGapThreshold = value;
}

unsigned Ptba::run() const noexcept(false) {
  auto result = result_from_run();

  /*
  {
    Ptba::dumpEigenGaps(result.eigenGaps, "eigengaps.txt");
    Ptba::dumpPerturbedEigenGaps(result.perturbedEigenGaps,
                                 "perturbed_eigengaps.txt");
  }
  */

  return result.significantIndices.size();
}

auto Ptba::result_from_run() const noexcept(false) -> PtbaResult {
  enum class PValueResult { significant, nonsignificant, inf, alternative };

  arma::mat dataEigenVecs;
  arma::vec dataEigenVals;
  arma::vec dataEigenGaps;
  arma::mat adjacency;
  auto initialData = *ringmapData;
  std::vector<unsigned> filteredToUnfilteredBases(
      initialData.getSequence().size());
  initialData.filterBases();

  {
    auto filteredData = initialData;
    filteredData.filterReads();
    for (auto [filteredIndex, unfilteredIndex] :
         filteredData.getFilteredToNonFilteredMap()) {
      assert(filteredIndex < filteredToUnfilteredBases.size());
      filteredToUnfilteredBases[filteredIndex] = unfilteredIndex;
    }
    if (filteredData.size() < minFilteredReads or
        filteredData.data().cols_size() < minBasesSize)
      return {};

    std::tie(dataEigenVecs, dataEigenVals, dataEigenGaps, adjacency) =
        calculateEigenGaps(filteredData);
    assert(dataEigenGaps.size() > 1);

    if (arma::all(dataEigenGaps == 0))
      return {};
  }

  assert(dataEigenGaps.size() > 1);
  assert(maxClusters > 0);
  const auto useful_eigengaps =
      std::min(static_cast<unsigned>(dataEigenGaps.size() / 2), maxClusters);

  if (useful_eigengaps == 0)
    return {};

#ifndef NDEBUG
  PerturbedEigengaps perturbed_eigengaps(useful_eigengaps);
  dataEigenGaps.resize(useful_eigengaps);
#else
  PerturbedEigengaps perturbed_eigengaps(dataEigenGaps.size());
#endif
  std::vector<WeibullFitter> weibull_fitters;
  weibull_fitters.reserve(dataEigenGaps.size());
  for (auto const &perturbed_eigengap : perturbed_eigengaps) {
    weibull_fitters.emplace_back(perturbed_eigengap);
  }
  std::vector<WeibullParams> weibull_params(dataEigenGaps.size(), {-1., -1.});
  std::vector<std::uint8_t> null_dist_is_valid(dataEigenGaps.size(), 0);

  unsigned eigenGapIndex = 0u;
  unsigned valid_eigengap_index = std::numeric_limits<unsigned>::max();
  assert(initialData.size() != 0);
  for (unsigned permutation = 0;
       permutation < maxPermutations and eigenGapIndex < useful_eigengaps and
       (valid_eigengap_index == std::numeric_limits<unsigned>::max() or
        eigenGapIndex <= valid_eigengap_index + extended_search_eigengaps);
       ++permutation) {
    RingmapData perturbedData = initialData;
    perturbedData.perturb();
    perturbedData.filterBases();
    perturbedData.filterReads();

    if (perturbedData.size() < minFilteredReads)
      return {};

    auto current_perturbed_eigengaps =
        std::get<2>(calculateEigenGaps(perturbedData));
    {
      auto perturbed_eigengaps_iter = std::begin(perturbed_eigengaps);
      auto const perturbed_eigengaps_end = std::end(perturbed_eigengaps);
      auto current_perturbed_eigengaps_iter =
          std::cbegin(current_perturbed_eigengaps);
      for (; perturbed_eigengaps_iter < perturbed_eigengaps_end;
           ++perturbed_eigengaps_iter, ++current_perturbed_eigengaps_iter) {
        perturbed_eigengaps_iter->emplace_back(
            std::max(*current_perturbed_eigengaps_iter, 1e-6));
      }
    }

    if (permutation < minPermutations)
      continue;

    assert(std::all_of(
        std::begin(perturbed_eigengaps), std::end(perturbed_eigengaps),
        [&](auto &&perturbed_eigengap) {
          return perturbed_eigengap.size() == perturbed_eigengaps[0].size();
        }));

    auto isPValueSignificative =
        [&, this](unsigned eigenGapIndex) -> std::pair<PValueResult, bool> {
      auto &&perturbed_eigengap =
          std::as_const(perturbed_eigengaps)[eigenGapIndex];

      if (not null_dist_is_valid[eigenGapIndex]) {
        double const eigengap_value = dataEigenGaps[eigenGapIndex];
        double const null_mean = std::accumulate(
            std::begin(perturbed_eigengap), std::end(perturbed_eigengap), 0.,
            [n = static_cast<double>(perturbed_eigengap.size())](
                double acc, double x) { return acc + x / n; });
        double const null_stddev = std::sqrt(std::accumulate(
            std::begin(perturbed_eigengap), std::end(perturbed_eigengap), 0.,
            [&, n = static_cast<double>(perturbed_eigengap.size() - 1)](
                double acc, double x) {
              return acc + std::pow(x - null_mean, 2.) / n;
            }));

        if (null_stddev < min_null_stddev) {
          bool const away_from_null_model = [&] {
            if (eigengap_value < null_mean)
              return eigengap_value < null_mean - null_stddev * 3.;
            else
              return eigengap_value > null_mean + null_stddev * 3.;
          }();

          if (away_from_null_model)
            return std::pair(PValueResult::significant, false);
          else if (eigengap_value >= null_mean - null_stddev * 3. and
                   eigengap_value <= null_mean + null_stddev * 3.)
            return std::pair(PValueResult::alternative, false);
          else
            return std::pair(PValueResult::nonsignificant, false);
        } else if (permutation >= alternative_check_permutations)
          return std::pair(PValueResult::significant, false);
      }

      auto &params = weibull_params[eigenGapIndex];
      params = weibull_fitters[eigenGapIndex].fit();
      boost::math::weibull_distribution<double> distribution(params.shape,
                                                             params.scale);

      if (not null_dist_is_valid[eigenGapIndex]) {
        null_dist_is_valid[eigenGapIndex] =
            is_distribution(distribution, perturbed_eigengap);
        if (not null_dist_is_valid[eigenGapIndex])
          return std::pair(PValueResult::inf, true);
      }

      double pValue;
      if (eigenGapIndex == 0)
        pValue = boost::math::cdf(distribution, dataEigenGaps[eigenGapIndex]);
      else {
        pValue = boost::math::cdf(boost::math::complement(
            distribution, dataEigenGaps[eigenGapIndex]));
      }

      if (pValue < alphaValue)
        return std::pair(PValueResult::significant, true);
      else if (pValue < betaValue)
        return std::pair(PValueResult::nonsignificant, true);
      else
        return std::pair(PValueResult::alternative, true);
    };

    if (eigenGapIndex == 0) {
      auto const [significativeness, evaluated_distribution] =
          isPValueSignificative(0);

      if (significativeness == PValueResult::inf) {
        if (permutation == maxPermutations - 1)
          return {0u,
                  std::move(dataEigenVecs),
                  std::move(dataEigenGaps),
                  std::move(perturbed_eigengaps),
                  {},
                  std::move(filteredToUnfilteredBases),
                  std::move(adjacency)};
        else
          continue;
      }

      if (not evaluated_distribution and
          significativeness != PValueResult::significant) {
        if (permutation < maxPermutations - 1)
          continue;
        else {
          auto const &perturbed_eigengap = perturbed_eigengaps[0];
          double shifted_mean;
          double first_stddev;
          {
            double const first_mean = std::accumulate(
                std::begin(perturbed_eigengap), std::end(perturbed_eigengap),
                0.,
                [n = static_cast<double>(perturbed_eigengap.size())](
                    double acc, double x) { return acc + x / n; });
            shifted_mean = first_mean * firstEigengapThreshold;

            first_stddev = std::sqrt(std::accumulate(
                std::begin(perturbed_eigengap), std::end(perturbed_eigengap),
                0.,
                [&, n = static_cast<double>(perturbed_eigengap.size() - 1)](
                    double acc, double x) {
                  return acc + std::pow(x - first_mean, 2.) / n;
                }));
          }

          if (significativeness == PValueResult::alternative and
              dataEigenGaps[0] < shifted_mean - first_stddev * 3.) {
            return {1u,
                    std::move(dataEigenVecs),
                    std::move(dataEigenGaps),
                    std::move(perturbed_eigengaps),
                    std::vector<unsigned>{0},
                    std::move(filteredToUnfilteredBases),
                    std::move(adjacency)};
          } else {
            return {0u,
                    std::move(dataEigenVecs),
                    std::move(dataEigenGaps),
                    std::move(perturbed_eigengaps),
                    {},
                    std::move(filteredToUnfilteredBases),
                    std::move(adjacency)};
          }
        }
      }

      if (evaluated_distribution) {
        double const first_mean =
            weibull_params[0].scale *
            std::tgamma(1. + 1. / weibull_params[0].shape);

        {
          PerturbedEigengap shiftedFirstEigengap(perturbed_eigengaps[0].size());
          std::transform(
              std::cbegin(perturbed_eigengaps[0]),
              std::cend(perturbed_eigengaps[0]),
              std::begin(shiftedFirstEigengap), [&](auto &&value) {
                return std::max(
                    value - (1. - firstEigengapThreshold) * first_mean, 1e-5);
              });
          assert(std::none_of(std::begin(shiftedFirstEigengap),
                              std::end(shiftedFirstEigengap), [](double value) {
                                return std::isinf(value) or std::isnan(value);
                              }));

          WeibullFitter fitter(shiftedFirstEigengap);
          auto shifted_params = fitter.fit();

          boost::math::weibull_distribution<double> shifted_distribution(
              shifted_params.shape, shifted_params.scale);

          if (double pValue =
                  boost::math::cdf(shifted_distribution, dataEigenGaps[0]);
              pValue > firstEigengapBetaValue) {
            return {0u,
                    std::move(dataEigenVecs),
                    std::move(dataEigenGaps),
                    std::move(perturbed_eigengaps),
                    {},
                    std::move(filteredToUnfilteredBases),
                    std::move(adjacency)};
          } else if (pValue >= alphaValue / 2. or
                     std::abs(first_mean - dataEigenGaps[0]) <
                         first_mean * (1. - firstEigengapThreshold)) {
            if (permutation < maxPermutations - 1)
              continue;
            else
              return {0u,
                      std::move(dataEigenVecs),
                      std::move(dataEigenGaps),
                      std::move(perturbed_eigengaps),
                      {},
                      std::move(filteredToUnfilteredBases),
                      std::move(adjacency)};
          }
        }
      } else {
        // Approximated evaluation
        auto const &perturbed_eigengap = perturbed_eigengaps[0];
        double shifted_mean;
        double first_stddev;
        {
          double const first_mean = std::accumulate(
              std::begin(perturbed_eigengap), std::end(perturbed_eigengap), 0.,
              [n = static_cast<double>(perturbed_eigengap.size())](
                  double acc, double x) { return acc + x / n; });
          shifted_mean = first_mean * firstEigengapThreshold;

          first_stddev = std::sqrt(std::accumulate(
              std::begin(perturbed_eigengap), std::end(perturbed_eigengap), 0.,
              [&, n = static_cast<double>(perturbed_eigengap.size() - 1)](
                  double acc, double x) {
                return acc + std::pow(x - first_mean, 2.) / n;
              }));
        }

        double const first_eigengap = dataEigenGaps[0];
        if (first_eigengap >= shifted_mean + first_stddev * 3.) {
          return {0u,
                  std::move(dataEigenVecs),
                  std::move(dataEigenGaps),
                  std::move(perturbed_eigengaps),
                  {},
                  std::move(filteredToUnfilteredBases),
                  std::move(adjacency)};
        } else if (first_eigengap >= shifted_mean - first_stddev * 3.) {
          if (permutation < maxPermutations - 1)
            continue;
          else
            return {0u,
                    std::move(dataEigenVecs),
                    std::move(dataEigenGaps),
                    std::move(perturbed_eigengaps),
                    {},
                    std::move(filteredToUnfilteredBases),
                    std::move(adjacency)};
        }
      }
    }

    auto const calc_eigengap_null_distance = [&](unsigned eigengap_index,
                                                 bool evaluated_distribution) {
      double mean, stddev;
      if (evaluated_distribution) {
        auto &params = weibull_params[eigengap_index];
        params = weibull_fitters[eigengap_index].fit();
        double const common = std::tgamma(1. + 1. / params.shape);
        mean = params.scale * common;
        stddev = params.scale * std::sqrt(std::tgamma(1. + 2. / params.shape) -
                                          std::pow(common, 2.));
      } else {
        auto const &perturbed_eigengap = perturbed_eigengaps[eigengap_index];
        mean = std::accumulate(
            std::begin(perturbed_eigengap), std::end(perturbed_eigengap), 0.,
            [n = static_cast<double>(perturbed_eigengap.size())](
                double acc, double x) { return acc + x / n; });
        stddev = std::sqrt(std::accumulate(
            std::begin(perturbed_eigengap), std::end(perturbed_eigengap), 0.,
            [&, n = static_cast<double>(perturbed_eigengap.size() - 1)](
                double acc, double x) {
              return acc + std::pow(x - mean, 2.) / n;
            }));
      }

      if (eigengap_index == 0)
        return std::max(mean - stddev - dataEigenGaps[eigengap_index], 0.);
      else
        return std::abs(dataEigenGaps[eigengap_index] - mean - stddev);
    };

    eigenGapIndex = std::max(eigenGapIndex, 1u);
    if (valid_eigengap_index == std::numeric_limits<unsigned>::max()) {
      valid_eigengap_index = 0u;
    }
    for (; eigenGapIndex < useful_eigengaps and
           eigenGapIndex <= valid_eigengap_index + extended_search_eigengaps;
         ++eigenGapIndex) {
      auto [result, evaluated_distribution] =
          isPValueSignificative(eigenGapIndex);

      if (result == PValueResult::inf or
          result == PValueResult::nonsignificant) {
        if (eigenGapIndex > valid_eigengap_index + extended_search_eigengaps or
            permutation < maxPermutations - 1)
          break;
      }

      if (result == PValueResult::significant) {
        const double diff =
            calc_eigengap_null_distance(eigenGapIndex, evaluated_distribution);
        if (diff >= eigenGapDiffAbsoluteThreshold) {
          double cum_diff = 0.;
          for (unsigned eigengap_index = 0; eigengap_index < eigenGapIndex;
               ++eigengap_index) {

            // This is a rough evaluation, we can use the sample mean and stddev
            cum_diff += calc_eigengap_null_distance(eigengap_index, false);
          }

          if (diff >= cum_diff * minEigenGapThreshold) {
            valid_eigengap_index = eigenGapIndex;
            continue;
          }
        }
      }
    }

    if (eigenGapIndex == useful_eigengaps or
        eigenGapIndex > valid_eigengap_index + extended_search_eigengaps) {
      assert(valid_eigengap_index != std::numeric_limits<unsigned>::max());
      std::vector<unsigned> significantIndices(valid_eigengap_index + 1);
      std::iota(std::begin(significantIndices), std::end(significantIndices),
                0u);

      return {
          valid_eigengap_index + 1,      std::move(dataEigenVecs),
          std::move(dataEigenGaps),      std::move(perturbed_eigengaps),
          std::move(significantIndices), std::move(filteredToUnfilteredBases),
          std::move(adjacency)};
    }
  }

  constexpr std::string_view eigenVecsFilename = "/tmp/eigenvecs_thrown.txt";
  constexpr std::string_view eigenGapsFilename = "/tmp/eigengaps_thrown.txt";
  constexpr std::string_view perturbedEigenGapsFilename =
      "/tmp/perturbed_eigengaps_thrown.txt";
  try {
    dumpEigenVecs(dataEigenVecs, eigenVecsFilename);
  } catch (const std::exception &) {
    std::cerr << "WARNING: eigenvecs data cannot be written\n";
  }

  try {
    dumpEigenGaps(dataEigenGaps, eigenGapsFilename);
  } catch (const std::exception &) {
    std::cerr << "WARNING: eigengaps data cannot be written\n";
  }

  try {
    dumpPerturbedEigenGaps(perturbed_eigengaps, perturbedEigenGapsFilename);
  } catch (const std::exception &) {
    std::cerr << "WARNING: perturbed eigengaps data cannot be written\n";
  }

  throw exception(std::string("cannot find the number of clusters. eigenvecs, "
                              "eigengaps and perturbed eigengaps written to ") +
                  std::string(eigenVecsFilename) + ", " +
                  std::string(eigenGapsFilename) + " and " +
                  std::string(perturbedEigenGapsFilename));
}

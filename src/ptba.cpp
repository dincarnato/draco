#include "ptba.hpp"
#include "draco.hpp"
#include "logger.hpp"
#include "results/transcript.hpp"
#include "weibull_fitter.hpp"

#include <algorithm>
#include <format>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <ranges>
#include <type_traits>
#include <variant>

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
      minEigenGapThreshold(args.min_eigengap_threshold()),
      eigenGapDiffAbsoluteThreshold(args.eigengap_diff_absolute_threshold()),
      alphaValue(args.alpha_value()), betaValue(args.beta_value()),
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
  PtbaResult result = run();
  dumpEigenVecs(result.eigenVecs, std::move(eigenVecsFilename));
  dumpEigenGaps(result.eigenGaps, std::move(eigenGapsFilename));
  dumpPerturbedEigenGaps(result.perturbedEigenGaps,
                         std::move(perturbedEigenGapsFilename));

  return result.nClusters;
}

std::pair<unsigned, std::vector<unsigned>>
Ptba::getNumberOfClustersAndSignificantIndices() const {
  PtbaResult result = run();
  return {result.nClusters, std::move(result.significantIndices)};
}

std::pair<unsigned, std::vector<unsigned>>
Ptba::getNumberOfClustersAndSignificantIndicesAndDumpData(
    std::string_view eigenGapsFilename,
    std::string_view perturbedEigenGapsFilename) const {
  PtbaResult result = run();
  dumpEigenGaps(result.eigenGaps, std::move(eigenGapsFilename));
  dumpPerturbedEigenGaps(result.perturbedEigenGaps,
                         std::move(perturbedEigenGapsFilename));

  return {result.nClusters, std::move(result.significantIndices)};
}

void Ptba::setMaxClusters(unsigned value) { maxClusters = value; }

std::vector<unsigned> Ptba::getAllSignificantEigenGapIndices() const {
  return run().significantIndices;
}

void Ptba::setMinEigenGapThreshold(double value) {
  minEigenGapThreshold = value;
}

auto Ptba::run() const noexcept(false) -> PtbaResult {
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

    if (filteredData.size() < minFilteredReads) {
      return {LogData{
          log_data::NotEnoughReads{.reads = filteredData.size(),
                                   .min_filtered_reads = minFilteredReads}}};
    } else if (filteredData.data().cols_size() < minBasesSize) {
      return {LogData{
          log_data::NotEnoughBases{.bases = filteredData.data().cols_size(),
                                   .min_bases = minBasesSize}}};
    }

    std::tie(dataEigenVecs, dataEigenVals, dataEigenGaps, adjacency) =
        calculateEigenGaps(filteredData);
    assert(dataEigenGaps.size() > 1);

    if (arma::all(dataEigenGaps == 0))
      return {LogData{log_data::AllZeroEigenGaps{}}};
  }

  assert(dataEigenGaps.size() > 1);
  assert(maxClusters > 0);
  const auto useful_eigengaps =
      std::min(static_cast<unsigned>(dataEigenGaps.size() / 2), maxClusters);

  if (useful_eigengaps == 0)
    return {LogData{log_data::NoUsefulEigenGaps{
        .eigengaps = dataEigenGaps.size(), .max_clusters = maxClusters}}};

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

  log_data::Permuting permutation_log;
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

    if (perturbedData.size() < minFilteredReads) {
      logger::on_debug_level([&]() {
        permutation_log.data.push_back(log_data::permuting::NotEnoughReads{
            .permutation = permutation,
            .filtered_reads = perturbedData.size(),
            .min_filtered_reads = minFilteredReads});
      });
      return {LogData{log_data::Permuting{std::move(permutation_log)}}};
    }

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

      // We have modified this part of the code to calcualte a two-tailed
      // p-value, after seeing cases of switches with the null model and the
      // eigengap inverted if (eigenGapIndex == 0)
      double cdfValue =
          boost::math::cdf(distribution, dataEigenGaps[eigenGapIndex]);
      double pValue = 2 * std::min(cdfValue, 1 - cdfValue);

      if (pValue < alphaValue)
        return std::pair(PValueResult::significant, true);
      else if (pValue < betaValue)
        return std::pair(PValueResult::nonsignificant, true);
      else
        return std::pair(PValueResult::alternative, true);
    };

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

      return std::max(std::abs(dataEigenGaps[eigengap_index] - mean) - stddev,
                      0.);
    };

    eigenGapIndex = std::max(eigenGapIndex, 1u);
    if (valid_eigengap_index == std::numeric_limits<unsigned>::max()) {
      valid_eigengap_index = 0u;
    }
    for (; eigenGapIndex < useful_eigengaps and
           eigenGapIndex <= valid_eigengap_index + extended_search_eigengaps;
         ++eigenGapIndex) {
      auto result_and_evaluated_distribution =
          isPValueSignificative(eigenGapIndex);
      auto result = result_and_evaluated_distribution.first;
      auto evaluated_distribution = result_and_evaluated_distribution.second;

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
            logger::on_debug_level([&]() {
              permutation_log.data.push_back(
                  log_data::permuting::NewValidEigengap{
                      .permutation = permutation,
                      .previous_eigengap_index = valid_eigengap_index,
                      .new_eigengap_index = eigenGapIndex,
                      .distribution_is_evaluated = evaluated_distribution,
                      .eigengap_difference = diff,
                      .cumulative_difference = cum_diff,
                      .min_eigengap_threshold = minEigenGapThreshold,
                  });
            });

            valid_eigengap_index = eigenGapIndex;
            continue;
          }
        } else {
          logger::on_debug_level([&]() {
            permutation_log.data.push_back(
                log_data::permuting::SignificantEigengapLowDifference{
                    .permutation = permutation,
                    .eigengap_difference = diff,
                    .eigengap_diff_absolute_threshold =
                        eigenGapDiffAbsoluteThreshold,
                    .distribution_is_evaluated = evaluated_distribution,
                    .valid_eigengap_index = valid_eigengap_index,
                    .current_eigengap_index = eigenGapIndex,
                });
          });
        }
      } else {
        logger::on_debug_level([&]() {
          permutation_log.data.push_back(
              log_data::permuting::NotSignificantEigengap{
                  .permutation = permutation,
                  .valid_eigengap_index = valid_eigengap_index,
                  .current_eigengap_index = eigenGapIndex,
              });
        });
      }
    }

    if (eigenGapIndex == useful_eigengaps or
        eigenGapIndex > valid_eigengap_index + extended_search_eigengaps) {
      assert(valid_eigengap_index != std::numeric_limits<unsigned>::max());

      logger::on_debug_level([&]() {
        auto solution = ([&]() {
          if (eigenGapIndex == useful_eigengaps) {
            return log_data::permuting::solution_found::Solution{
                log_data::permuting::solution_found::UsefulEigengapIndex{}};
          } else {
            return log_data::permuting::solution_found::Solution{
                log_data::permuting::solution_found::OverExtendedSearch{
                    .valid_eigengap_index = valid_eigengap_index,
                    .extended_search_eigengaps = extended_search_eigengaps,
                }};
          }
        })();

        permutation_log.data.push_back(log_data::permuting::SolutionFound{
            .permutation = permutation,
            .eigengap_index = eigenGapIndex,
            .useful_eigengaps = useful_eigengaps,
            .solution = solution,
        });
      });

      std::vector<unsigned> significantIndices(valid_eigengap_index + 1);
      std::iota(std::begin(significantIndices), std::end(significantIndices),
                0u);

      return {
          valid_eigengap_index + 1,      std::move(dataEigenVecs),
          std::move(dataEigenGaps),      std::move(perturbed_eigengaps),
          std::move(significantIndices), std::move(filteredToUnfilteredBases),
          std::move(adjacency),          std::move(permutation_log)};
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

void print_log_data(LogData const &log_data, Window const &window,
                    size_t window_index, unsigned window_size,
                    results::Transcript const &transcript) noexcept {
  std::visit(
      [&](auto const &log_data) {
        using variant = std::decay_t<decltype(log_data)>;
        if constexpr (std::is_same_v<variant, log_data::NotEnoughReads>) {
          logger::debug(
              "Transcript {}, window {} (bases {}-{}) is skipped by "
              "initial filtering: reads ({}) < min_filtered_reads ({}).",
              transcript.name, window_index + 1, window.start_base + 1,
              window.start_base + window_size, log_data.reads,
              log_data.min_filtered_reads);
        } else if constexpr (std::is_same_v<variant,
                                            log_data::NotEnoughBases>) {
          logger::debug("Transcript {}, window {} (bases {}-{}) is skipped by "
                        "initial filtering: bases ({}) < min_bases ({}).",
                        transcript.name, window_index + 1,
                        window.start_base + 1, window.start_base + window_size,
                        log_data.bases, log_data.min_bases);
        } else if constexpr (std::is_same_v<variant,
                                            log_data::AllZeroEigenGaps>) {
          logger::debug("Transcript {}, window {} (bases {}-{}) is skipped "
                        "because all eigengaps are zero.",
                        transcript.name, window_index + 1,
                        window.start_base + 1, window.start_base + window_size);
        } else if constexpr (std::is_same_v<variant,
                                            log_data::NoUsefulEigenGaps>) {
          logger::debug("Transcript {}, window {} (bases {}-{}) is skipped "
                        "because the number of eigengaps ({}) or the max "
                        "number of clusters ({}) is less than 1",
                        transcript.name, window_index + 1,
                        window.start_base + 1, window.start_base + window_size,
                        log_data.eigengaps, log_data.max_clusters);
        } else if constexpr (std::is_same_v<variant, log_data::Permuting>) {
          std::string message = std::format(
              "Transcript {}, window {} (bases {}-{}):", transcript.name,
              window_index + 1, window.start_base + 1,
              window.start_base + window_size);

          if (log_data.data.empty()) {
            logger::debug("{} finished PTBA without any information", message);
            return;
          }

          auto print_permutation_log_data = [&](auto const
                                                    &permutation_log_data) {
            std::visit(
                [&](auto const &permutation_log_data) {
                  using variant = std::decay_t<decltype(permutation_log_data)>;
                  if constexpr (std::is_same_v<
                                    variant,
                                    log_data::permuting::NotEnoughReads>) {
                    std::format_to(std::back_inserter(message),
                                   " does not have enough reads ({}, min = "
                                   "{}) on permutation {}",
                                   permutation_log_data.filtered_reads,
                                   permutation_log_data.min_filtered_reads,
                                   permutation_log_data.permutation + 1);
                  } else if constexpr (std::is_same_v<variant,
                                                      log_data::permuting::
                                                          NewValidEigengap>) {
                    std::format_to(
                        std::back_inserter(message),
                        " new valid eigengap found (index = {}, previous "
                        "index = {}) on permutation {}, eigengap "
                        "difference ({}) >= cumulative difference ({}) * "
                        "min eigengap threshold ({})",
                        permutation_log_data.new_eigengap_index,
                        permutation_log_data.previous_eigengap_index,
                        permutation_log_data.permutation + 1,
                        permutation_log_data.eigengap_difference,
                        permutation_log_data.cumulative_difference,
                        permutation_log_data.min_eigengap_threshold);
                  } else if constexpr (
                      std::is_same_v<variant,
                                     log_data::permuting::
                                         SignificantEigengapLowDifference>) {
                    std::format_to(
                        std::back_inserter(message),
                        " eigengap with index {} is significant on permutation "
                        "{}, but the absolute eigengap difference ({}) is too "
                        "low (threshold = {}) -- valid eigengap index = {}, "
                        "distribution is {}evaluated",
                        permutation_log_data.current_eigengap_index,
                        permutation_log_data.permutation + 1,
                        permutation_log_data.eigengap_difference,
                        permutation_log_data.eigengap_diff_absolute_threshold,
                        permutation_log_data.valid_eigengap_index,
                        permutation_log_data.distribution_is_evaluated
                            ? ""
                            : "not ");
                  } else if constexpr (std::is_same_v<
                                           variant,
                                           log_data::permuting::
                                               NotSignificantEigengap>) {
                    std::format_to(
                        std::back_inserter(message),
                        " eigengap with index {} is not significant on "
                        "permutation {} -- valid eigengap index = {}",
                        permutation_log_data.current_eigengap_index,
                        permutation_log_data.permutation + 1,
                        permutation_log_data.valid_eigengap_index);
                  } else if constexpr (std::is_same_v<variant,
                                                      log_data::permuting::
                                                          SolutionFound>) {
                    std::visit(
                        [&](auto const &solution) {
                          using variant = std::decay_t<decltype(solution)>;
                          if constexpr (
                              std::is_same_v<
                                  variant, log_data::permuting::solution_found::
                                               UsefulEigengapIndex>) {
                            std::format_to(
                                std::back_inserter(message),
                                " found a solution in permutation {}, with "
                                "eigengap index ({}) == useful eigengaps ({}) ",
                                permutation_log_data.permutation + 1,
                                permutation_log_data.eigengap_index,
                                permutation_log_data.useful_eigengaps);

                          } else if constexpr (
                              std::is_same_v<
                                  variant, log_data::permuting::solution_found::
                                               OverExtendedSearch>) {
                            std::format_to(
                                std::back_inserter(message),
                                " found a solution in permutation {}, with "
                                "eigengap index {}, reached checking eigengap "
                                "index {} (> valid index + {}) without any "
                                "useful "
                                "information",
                                permutation_log_data.permutation + 1,
                                solution.valid_eigengap_index,
                                permutation_log_data.eigengap_index,
                                solution.extended_search_eigengaps);

                          } else {
                            static_assert(false, "non-exhaustive visitor");
                          }
                        },
                        permutation_log_data.solution);
                  } else {
                    static_assert(false, "non-exhaustive visitor");
                  }
                },
                permutation_log_data);
          };

          print_permutation_log_data(log_data.data[0]);
          std::ranges::for_each(log_data.data | std::views::drop(1),
                                [&](auto const &permutation_log_data) {
                                  message.push_back(';');
                                  print_permutation_log_data(
                                      permutation_log_data);
                                });

          logger::debug("{}", message);
        } else {
          static_assert(false, "non-exhaustive visitor");
        }
      },
      log_data);
}

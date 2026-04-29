#include "reassignment.hpp"
#include "compact_ringmap.hpp"
#include "draco.hpp"
#include "expectation_maximization.hpp"
#include "logger.hpp"
#include "read_clusters_assignments.hpp"
#include "results/transcript.hpp"
#include "results/window.hpp"
#include "ringmap_data.hpp"

#include <algorithm>
#include <format>
#include <functional>
#include <iterator>
#include <random>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>

void Reassignment::reassign_reads_with_weights(
    bool allow_empty_patterns) const {
  for (auto &filtered_ringmap : filtered_ringmaps) {
    filtered_ringmap.filterReads();
  }

  std::ranges::for_each(
      std::views::zip(std::views::iota(0uz), filtered_ringmaps,
                      replicates_splitted_ringmaps,
                      transcript_result->windows |
                          std::views::transform(
                              [&](auto &replicate_windows) -> decltype(auto) {
                                return (*replicate_windows)[window_index];
                              })) |
          std::views::take_while([&](auto const &) { return *stop; }),
      [&](auto &&tuple) {
        auto replicate_index = std::get<0>(tuple);
        auto const &filtered_ringmap = std::get<1>(tuple);
        auto const &ringmap = std::get<2>(tuple)[window_index];
        auto &window = std::get<3>(tuple);

        auto &&fractions_result = filtered_ringmap.fractionReadsByWeights(
            window.weighted_clusters, window_size,
            args->skip_ambiguous_assignments());

        reassignment::HandleFractionedReads{
            .window = &window,
            .filtered_ringmap = &filtered_ringmap,
            .ringmap = &ringmap,
            .stop = stop,
            .transcript_result = transcript_result,
            .args = args,
            .replicate_index = replicate_index,
            .window_index = window_index,
            .window_size = window_size,
            .ptba_on_replicate_results = ptba_on_replicate_results,
            .windows_max_clusters_constraints =
                windows_max_clusters_constraints,
            .fractions_result = std::move(fractions_result),
            .allow_empty_patterns = allow_empty_patterns,
        }();
      });
}

void Reassignment::reweight_and_reassign_with_expectation_maximization() const {
  std::mt19937 rng(std::random_device{}());
  std::vector<std::uint32_t> assignments_per_cluster;
  std::vector<double> buffer;
  std::vector<std::uint32_t> mapped_rows;
  std::vector<std::uint32_t> clusters_reads_count;

  std::ranges::for_each(
      std::views::zip(std::views::iota(0uz), replicates_splitted_ringmaps,
                      transcript_result->windows |
                          std::views::transform(
                              [&](auto &replicate_windows) -> decltype(auto) {
                                return (*replicate_windows)[window_index];
                              })) |
          std::views::take_while([&](auto const &) { return *stop; }),
      [&](auto &&tuple) {
        auto replicate_index = std::get<0>(tuple);
        auto const &ringmap = std::get<1>(tuple)[window_index];
        auto &window = std::get<2>(tuple);

        auto fractions_result =
            reweight_and_reassign_with_expectation_maximization_iteration(
                reassignment::
                    ReweightAndReassignWithExpectationMaximizationIteration{
                        .ringmap = &ringmap,
                        .window = &window,
                        .rng = &rng,
                        .assignments_per_cluster = &assignments_per_cluster,
                        .buffer = &buffer,
                        .mapped_rows = &mapped_rows,
                        .clusters_reads_count = &clusters_reads_count,
                        .replicate_index = replicate_index,
                    });

        reassignment::HandleFractionedReads{
            .window = &window,
            .filtered_ringmap = nullptr,
            .ringmap = &ringmap,
            .stop = stop,
            .transcript_result = transcript_result,
            .args = this->args,
            .replicate_index = replicate_index,
            .window_index = window_index,
            .window_size = window_size,
            .ptba_on_replicate_results = ptba_on_replicate_results,
            .windows_max_clusters_constraints =
                windows_max_clusters_constraints,
            .fractions_result = std::move(fractions_result),
            .allow_empty_patterns = true,
        }();
      });
}

reassignment::FractionResult
Reassignment::reweight_and_reassign_with_expectation_maximization_iteration(
    reassignment::ReweightAndReassignWithExpectationMaximizationIteration args)
    const {
  auto &weighted_clusters = args.window->weighted_clusters;
  CompactRingmap compact_ringmap(args.ringmap->data());
  ExpectationMaximization expectation_maximization(
      compact_ringmap, weighted_clusters, *this->args, *args.rng);
  auto em_result = expectation_maximization.run();
  std::visit(
      [&](auto &&convergence) {
        using T = std::remove_cvref_t<decltype(convergence)>;
        if constexpr (std::is_same_v<T, expectation_maximization::Converged>) {
          logger::debug("Expectation-maximization on transcript {}, "
                        "replicate {}, window {} (bases {}-{}) converged "
                        "to log-likelihood of {} after {} iterations",
                        transcript_result->name, args.replicate_index + 1,
                        window_index + 1, args.window->begin_index + 1,
                        args.window->end_index, em_result.log_likelihood,
                        convergence.after_iterations);
        } else if constexpr (std::is_same_v<
                                 T, expectation_maximization::MaxIterations>) {
          logger::warn("Expectation-maximization on transcript {}, "
                       "replicate {}, window {} (bases {}-{}) reached max "
                       "iterations without converging (log-likelihood = {})",
                       transcript_result->name, args.replicate_index + 1,
                       window_index + 1, args.window->begin_index + 1,
                       args.window->end_index, em_result.log_likelihood);
        } else {
          static_assert(false, "unreachable");
        }
      },
      em_result.convergence);

  RingmapData::clusters_assignment_type clusters_assignment;
  RingmapData::clusters_pattern_type patterns(
      weighted_clusters.getClustersSize(),
      RingmapData::cluster_pattern_type(weighted_clusters.getElementsSize(),
                                        0));
  args.clusters_reads_count->resize(weighted_clusters.getClustersSize());
  std::ranges::fill(*args.clusters_reads_count, static_cast<std::uint32_t>(0));

  args.assignments_per_cluster->resize(weighted_clusters.getClustersSize());
  args.buffer->resize(weighted_clusters.getClustersSize());

#ifndef NDEBUG
  std::uint32_t skipped_reads = 0;
#endif

  auto rows_iter = std::ranges::begin(compact_ringmap);
  while (rows_iter != std::ranges::end(compact_ringmap)) {
    auto rows_range_end = std::ranges::find_if(
        std::ranges::next(rows_iter), std::ranges::end(compact_ringmap),
        [&](auto &&row) {
          return not std::ranges::equal(row.indices(), (*rows_iter).indices());
        });

    CompactRingmapRange rows(rows_iter, rows_range_end);
    expectation_maximization.read_assignment(
        rows, *args.assignments_per_cluster, *args.buffer, *args.rng);
    auto mapped_rows_size =
        std::ranges::fold_left(rows | std::views::transform([](auto &&row) {
                                 return std::size(row.mapped_rows());
                               }),
                               0uz, std::plus{});

    args.mapped_rows->resize(mapped_rows_size);

    auto args_mapped_rows_iter = std::ranges::begin(*args.mapped_rows);
    for (auto &&row : rows) {
      args_mapped_rows_iter =
          std::ranges::copy(row.mapped_rows(), args_mapped_rows_iter).out;
    }
    std::ranges::shuffle(*args.mapped_rows, *args.rng);
    rows_iter = rows_range_end;

    std::uint32_t used_rows = 0;
    for (auto [cluster_index, assignments_count, cluster_patterns] :
         std::views::zip(std::views::iota(static_cast<std::uint8_t>(0)),
                         *args.assignments_per_cluster, patterns)) {
      auto cluster_mapped_rows =
          std::span(std::next(std::ranges::begin(*args.mapped_rows), used_rows),
                    assignments_count);
      std::ranges::sort(cluster_mapped_rows);

      for (auto read_index : cluster_mapped_rows) {
        auto read = args.ringmap->data().row(read_index);
        if (read.original_end_index() - read.original_begin_index() <
            window_size) {
#ifndef NDEBUG
          ++skipped_reads;
#endif
          logger::trace("Skipping read {}-{}, shorter than a window size ",
                        read.original_begin_index() + 1,
                        read.original_end_index());
          continue;
        }

        auto cluster_assignment_iter = ([&] {
          if (auto iter = clusters_assignment.find(read.modifiedIndices());
              iter != std::ranges::end(clusters_assignment)) {
            return iter;
          } else {
            return clusters_assignment
                .emplace(read.modifiedIndices(),
                         ReadClustersAssignments(
                             weighted_clusters.getClustersSize()))
                .first;
          }
        })();

        cluster_assignment_iter->second.cluster(cluster_index)
            .push_back(read_index);

        for (auto modified_index : read.modifiedIndices()) {
          cluster_patterns[modified_index] += 1;
        }
        (*args.clusters_reads_count)[cluster_index] += 1;
      }

      used_rows += assignments_count;
    }
  }
  // Check if we mapped in clusters_assignments all the reads compacted
  // in compact_ringmap
  assert(std::ranges::fold_left(
             clusters_assignment | std::views::transform([](auto &&tuple) {
               auto const &read_clusters_assignments = std::get<1>(tuple);
               return std::ranges::fold_left(
                   read_clusters_assignments.clusters() |
                       std::views::transform(
                           [](auto &&cluster) { return std::size(cluster); }),
                   static_cast<std::uint32_t>(0), std::plus{});
             }),
             static_cast<std::uint32_t>(0), std::plus{}) +
             skipped_reads ==
         std::ranges::fold_left(
             compact_ringmap |
                 std::views::transform([](auto &&row) { return row.count(); }),
             static_cast<std::uint32_t>(0), std::plus{}));

  RingmapData::clusters_fraction_type fractions(
      weighted_clusters.getClustersSize());
  auto total_reads_count = std::ranges::fold_left(
      *args.clusters_reads_count, static_cast<std::uint32_t>(0), std::plus{});
  std::ranges::transform(*args.clusters_reads_count,
                         std::ranges::begin(fractions),
                         [&](auto cluster_reads_count) {
                           return static_cast<double>(cluster_reads_count) /
                                  static_cast<double>(total_reads_count);
                         });

  return std::tuple{std::move(fractions), std::move(patterns),
                    std::move(clusters_assignment)};
}

namespace reassignment {

void HandleFractionedReads::operator()() {
  std::tie(window->fractions, window->patterns, std::ignore) =
      std::move(fractions_result);
  assert(window->fractions.size() > 1 or window->fractions.empty() or
         window->fractions[0] >= 0.01);

  bool const redundand_patterns = [&] {
    if (!allow_empty_patterns and
        std::ranges::any_of(*window->patterns, [](auto const &pattern) {
          return std::ranges::all_of(pattern,
                                     [](auto value) { return value == 0; });
        })) {
      throw std::runtime_error(std::format(
          "A pattern for window {}-{} of transcript {} contains only zeros. "
          "This should never happen.",
          window->begin_index, window->end_index, transcript_result->name));
    }

    auto patterns_iter = std::cbegin(*window->patterns);
    auto const patterns_end = std::cend(*window->patterns);

    for (; patterns_iter < patterns_end; ++patterns_iter) {
      auto &&cur_pattern = *patterns_iter;
      auto const begin_cur_pattern = std::cbegin(cur_pattern);
      auto const end_cur_pattern = std::cend(cur_pattern);

      if (std::any_of(std::next(patterns_iter), patterns_end,
                      [&](auto &&next_pattern) {
                        return std::equal(begin_cur_pattern, end_cur_pattern,
                                          std::cbegin(next_pattern),
                                          std::cend(next_pattern));
                      })) {
        return true;
      }
    }

    return false;
  }();

  if (redundand_patterns or
      std::ranges::any_of(
          window->fractions,
          [min_cluster_fraction = args->minimum_cluster_fraction()](
              auto &&fraction) { return fraction < min_cluster_fraction; })) {
    *stop = false;
    auto const result_window_begin = window->begin_index;
    auto const result_window_end = window->end_index;

    auto cause = ([&] {
      if (redundand_patterns) {
        return "a redundant weights pattern is found";
      } else {
        return "at least one fraction is below the minimum "
               "threshold";
      }
    })();
    logger::debug("Transcript {}, replicate {}, window {} (bases {}-{}), "
                  "reducing number of clusters from {} to {} because {}",
                  transcript_result->name, replicate_index + 1,
                  window_index + 1, result_window_begin + 1, result_window_end,
                  window->fractions.size(), window->fractions.size() - 1,
                  cause);
    assert(window->fractions.size() > 1);
    auto const new_clusters_constraint =
        static_cast<unsigned>(window->fractions.size() - 1);

    for (auto &replicate_result : ptba_on_replicate_results) {
      for (auto &&[window, window_constraint] : std::views::zip(
               replicate_result.windows, windows_max_clusters_constraints)) {
        if (window.start_base >= result_window_begin and
            window.start_base + window_size <= result_window_end) {

          window_constraint = new_clusters_constraint;
        }
      }
    }
  }

  if (not *stop)
    return;

  if (filtered_ringmap != nullptr) {
    *window->patterns = filtered_ringmap->remapPatterns(*window->patterns);
  }

  assign_reads_to_clusters(*window, std::move(std::get<2>(fractions_result)),
                           *ringmap, filtered_ringmap);

  if (not args->assignments_dump_directory().empty()) {
    std::optional<std::size_t> usable_replicate_index;
    if (std::size(ptba_on_replicate_results) > 1) {
      usable_replicate_index.emplace(replicate_index);
    }
    dump_assignments(*transcript_result, *window, *ringmap,
                     usable_replicate_index,
                     args->assignments_dump_directory());
  }

  if (std::all_of(std::cbegin(*window->patterns), std::cend(*window->patterns),
                  [](auto &&pattern) {
                    return std::all_of(std::cbegin(pattern), std::cend(pattern),
                                       [](auto &&value) { return value == 0; });
                  })) {
    window->patterns = std::nullopt;
    window->bases_coverages = std::nullopt;
  }
}

} // namespace reassignment

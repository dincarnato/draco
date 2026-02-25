#include "reassignment.hpp"
#include "draco.hpp"
#include "logger.hpp"
#include "results/transcript.hpp"
#include "ringmap_data.hpp"

#include <algorithm>
#include <format>
#include <iterator>
#include <ranges>
#include <tuple>
#include <utility>

void Reassignment::reassign_reads_with_weights() const {
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
        std::tie(window.fractions, window.patterns, std::ignore) =
            std::move(fractions_result);
        assert(window.fractions.size() > 1 or window.fractions.empty() or
               window.fractions[0] >= 0.01);

        bool const redundand_patterns = [&] {
          if (!allow_empty_patterns and
              std::ranges::any_of(*window.patterns, [](auto const &pattern) {
                return std::ranges::all_of(
                    pattern, [](auto value) { return value == 0; });
              })) {
            throw std::runtime_error(std::format(
                "A pattern for window {}-{} of "
                "transcript {} contains "
                "only zeros. This should never happen.",
                window.begin_index, window.end_index, transcript_result->name));
          }

          auto patterns_iter = std::cbegin(*window.patterns);
          auto const patterns_end = std::cend(*window.patterns);

          for (; patterns_iter < patterns_end; ++patterns_iter) {
            auto &&cur_pattern = *patterns_iter;
            auto const begin_cur_pattern = std::cbegin(cur_pattern);
            auto const end_cur_pattern = std::cend(cur_pattern);

            if (std::any_of(std::next(patterns_iter), patterns_end,
                            [&](auto &&next_pattern) {
                              return std::equal(begin_cur_pattern,
                                                end_cur_pattern,
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
                window.fractions,
                [min_cluster_fraction =
                     args->minimum_cluster_fraction()](auto &&fraction) {
                  return fraction < min_cluster_fraction;
                })) {
          *stop = false;
          auto const result_window_begin = window.begin_index;
          auto const result_window_end = window.end_index;

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
                        window_index + 1, result_window_begin + 1,
                        result_window_end, window.fractions.size(),
                        window.fractions.size() - 1, cause);
          assert(window.fractions.size() > 1);
          auto const new_clusters_constraint =
              static_cast<unsigned>(window.fractions.size() - 1);

          for (auto &replicate_result : ptba_on_replicate_results) {
            for (auto &&[window, window_constraint] :
                 std::views::zip(replicate_result.windows,
                                 windows_max_clusters_constraints)) {
              if (window.start_base >= result_window_begin and
                  window.start_base + window_size <= result_window_end) {

                window_constraint = new_clusters_constraint;
              }
            }
          }
        }

        if (not *stop)
          return;

        *window.patterns = filtered_ringmap.remapPatterns(*window.patterns);

        assign_reads_to_clusters(window,
                                 std::move(std::get<2>(fractions_result)),
                                 ringmap, filtered_ringmap);

        if (not args->assignments_dump_directory().empty()) {
          std::optional<std::size_t> usable_replicate_index;
          if (std::size(ptba_on_replicate_results) > 1) {
            usable_replicate_index.emplace(replicate_index);
          }
          dump_assignments(*transcript_result, window, ringmap,
                           usable_replicate_index,
                           args->assignments_dump_directory());
        }

        if (std::all_of(std::cbegin(*window.patterns),
                        std::cend(*window.patterns), [](auto &&pattern) {
                          return std::all_of(
                              std::cbegin(pattern), std::cend(pattern),
                              [](auto &&value) { return value == 0; });
                        })) {
          window.patterns = std::nullopt;
          window.bases_coverages = std::nullopt;
        }
      });
}


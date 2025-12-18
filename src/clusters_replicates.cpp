#include "clusters_replicates.hpp"
#include "fmt/base.h"
#include "logger.hpp"
#include "results/transcript.hpp"
#include "weighted_clusters.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

namespace clusters_replicates {
struct PermutationsFormatter {
  ReplicatesClustersPermutations const *replicates_clusters_permutations;
};
} // namespace clusters_replicates

template <> struct fmt::formatter<clusters_replicates::PermutationsFormatter> {
  fmt::format_context::iterator
  format(const clusters_replicates::PermutationsFormatter &f,
         fmt::format_context &ctx) const {
    auto out_iter = ctx.out();

    auto format_cluster_replicates = [&](auto &&cluster_replicates) {
      auto first_replicate_iter = std::ranges::begin(cluster_replicates);
      if (first_replicate_iter != std::ranges::end(cluster_replicates)) {
        ctx.advance_to(fmt::format_to(out_iter, "{}", *first_replicate_iter));

        for (auto &&cluster_replicate :
             cluster_replicates | std::views::drop(1)) {
          ctx.advance_to(fmt::format_to(out_iter, "-{}", cluster_replicate));
        }
      }
    };

    ctx.advance_to(fmt::format_to(out_iter, "["));
    auto &&clusters_replicates = f.replicates_clusters_permutations->clusters();
    if (not std::empty(clusters_replicates)) {
      if (auto first_cluster_replicates_iter =
              std::ranges::begin(clusters_replicates);
          first_cluster_replicates_iter !=
          std::ranges::end(clusters_replicates)) {
        format_cluster_replicates(*first_cluster_replicates_iter);

        for (auto &&cluster_replicates :
             clusters_replicates | std::views::drop(1)) {
          ctx.advance_to(fmt::format_to(out_iter, ", "));
          format_cluster_replicates(cluster_replicates);
        }
      }
    }

    return fmt::format_to(out_iter, "]");
  }

  constexpr auto parse(fmt::format_parse_context &ctx) { return ctx.begin(); }
};

namespace clusters_replicates {
void reorder_best_permutation(
    std::vector<WeightedClusters> &replicates_clusters,
    results::Transcript const &transcript, unsigned window_index,
    double distance_warning_threshold) {
  std::size_t const n_replicates = std::size(replicates_clusters);
  if (n_replicates <= 1) {
    return;
  }

  logger::trace("Searching the best permutation of clusters across replicates "
                "for transcript {} on window with index {}",
                transcript.name, window_index);

  PermutationsDistances permutation_distances(replicates_clusters);
  ReplicatesClustersPermutations replicates_clusters_permutations(
      permutation_distances.replicates(), permutation_distances.clusters());

  auto normalize_distance = [&](double distance) {
    auto n_bases = replicates_clusters[0].getElementsSize();
    return distance / (static_cast<double>(n_bases) *
                       static_cast<double>(
                           permutation_distances.replicates_combinationns()) *
                       static_cast<double>(permutation_distances.clusters()));
  };

  double best_clusters_distance = std::numeric_limits<double>::infinity();
  auto best_replicates_clusters_permutations = replicates_clusters_permutations;
  for (;;) {
    // TODO: maybe we should cache part of the calculation
    double total_clusters_distance = 0.;
    for (std::uint8_t replicate_1 = 0;
         replicate_1 < permutation_distances.replicates(); ++replicate_1) {
      auto replicate_1_clusters_permutations =
          replicates_clusters_permutations[replicate_1];

      for (std::uint8_t replicate_2 = replicate_1 + 1;
           replicate_2 < permutation_distances.replicates(); ++replicate_2) {
        auto replicate_2_clusters_permutations =
            replicates_clusters_permutations[replicate_2];
        auto replicates_distances =
            permutation_distances.replicates_pair_distances(replicate_1,
                                                            replicate_2);

        auto current_clusters_distance = replicates_distances.clusters_distance(
            replicate_1_clusters_permutations,
            replicate_2_clusters_permutations);
        total_clusters_distance += current_clusters_distance;
      }
    }

    if (total_clusters_distance < best_clusters_distance) {
      logger::on_trace_level([&] {
        logger::trace("New best permutation for transcript {} on window with "
                      "index {}; old distance was {} (normalized {}), new one "
                      "is {} (normalized {}), permutation indices are {}",
                      transcript.name, window_index, best_clusters_distance,
                      normalize_distance(best_clusters_distance),
                      total_clusters_distance,
                      normalize_distance(total_clusters_distance),
                      PermutationsFormatter{&replicates_clusters_permutations});
      });

      best_clusters_distance = total_clusters_distance;
      std::ranges::copy(
          replicates_clusters_permutations.indices(),
          std::ranges::begin(best_replicates_clusters_permutations.indices()));
    }

    std::uint8_t replicate_index = permutation_distances.replicates() - 1;
    for (;;) {
      if (std::ranges::next_permutation(
              replicates_clusters_permutations[replicate_index])
              .found) {
        break;
      }

      // We don't need to create different combinations for the first replicate
      if (replicate_index <= 1) {
        goto finished;
      } else {
        --replicate_index;
      }
    }
  }
finished:

  logger::on_debug_level([&] {
    logger::debug(
        "Best permutation for transcript {} on window with "
        "index {} has distance {} (normalized {}), permutations indices are {}",
        transcript.name, window_index, best_clusters_distance,
        normalize_distance(best_clusters_distance),
        PermutationsFormatter{&best_replicates_clusters_permutations});
  });

  assert(not std::isinf(best_clusters_distance));
  logger::on_warn_level([&] {
    auto normalized_distance = normalize_distance(best_clusters_distance);
    if (normalized_distance >= distance_warning_threshold) {
      logger::warn(
          "Best permutation for transcript {} on window with "
          "index {} has normalized distance {} (non-normalized {}) that is "
          "greater or equal than the warning threshold level ({}), "
          "permutations indices are {}",
          transcript.name, window_index, normalized_distance,
          best_clusters_distance, distance_warning_threshold,
          PermutationsFormatter{&best_replicates_clusters_permutations});
    }
  });

  std::ranges::for_each(
      std::views::zip(best_replicates_clusters_permutations,
                      replicates_clusters) |
          std::views::drop(1),
      [&](auto &&tuple) {
        auto &&[replicate_permutations, replicate_weighted_clusters] = tuple;

        for (std::uint8_t cluster_index = 0;
             cluster_index < permutation_distances.clusters() - 1;
             ++cluster_index) {

          auto &current_cluster_index = replicate_permutations[cluster_index];
          if (cluster_index != current_cluster_index) {
            auto swapping_cluster_iter = std::ranges::find(
                replicate_permutations | std::views::drop(cluster_index + 1),
                cluster_index);
            assert(swapping_cluster_iter !=
                   std::ranges::end(replicate_permutations));
            replicate_weighted_clusters.swap_clusters(cluster_index,
                                                      current_cluster_index);
            std::swap(current_cluster_index, *swapping_cluster_iter);
          }
        }
        assert(std::ranges::is_sorted(replicate_permutations));
      });
}

PermutationsDistances::PermutationsDistances(
    std::vector<WeightedClusters> &replicates_clusters)
    : replicates_(([&] {
        auto replicates = std::size(replicates_clusters);
        if (replicates > std::numeric_limits<std::uint8_t>::max()) {
          throw std::runtime_error("number of replicates is too high");
        }
        return static_cast<std::uint8_t>(replicates);
      })()),
      clusters_(([&] {
        if (replicates_clusters.empty()) {
          return static_cast<std::uint8_t>(0);
        }

        auto clusters = replicates_clusters[0].getClustersSize();
        if (clusters > std::numeric_limits<std::uint8_t>::max()) {
          throw std::runtime_error("number of clusters is too high");
        }
        return static_cast<std::uint8_t>(clusters);
      })()),
      replicates_combinations_(
          static_cast<std::uint16_t>(replicates_ * (replicates_ - 1) / 2)),
      distances_(distances_size(replicates_, clusters_)) {
  std::size_t distance_index = 0;
  for (std::uint8_t replicate_a_index = 0; replicate_a_index < replicates_;
       ++replicate_a_index) {
    auto const &replicate_a = replicates_clusters[replicate_a_index];

    for (std::uint8_t replicate_b_index = replicate_a_index + 1;
         replicate_b_index < replicates_; ++replicate_b_index) {
      auto const &replicate_b = replicates_clusters[replicate_b_index];

      for (auto replicate_a_clusters = replicate_a.clusters();
           auto &&cluster_a : replicate_a_clusters) {
        for (auto replicate_b_clusters = replicate_b.clusters();
             auto &&cluster_b : replicate_b_clusters) {
          auto distance = std::ranges::fold_left(
              std::views::zip(cluster_a, cluster_b) |
                  std::views::transform([](auto &&weights) {
                    return std::abs(std::get<0>(weights) -
                                    std::get<1>(weights));
                  }),
              0., std::plus{});
          distances_[distance_index++] = distance;
        }
      }
    }
  }
}

} // namespace clusters_replicates

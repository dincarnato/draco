#pragma once

#include "clusters_traits.hpp"
#include "graph_cut.hpp"
#include "logger.hpp"

#include <algorithm>
#include <armadillo>
#include <array>
#include <cassert>
#include <concepts>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <ranges>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

template <typename Fun, typename Ret = decltype(std::declval<Fun>()(
                            std::declval<const arma::mat &>()))>
static constexpr void checkGraphFunCallable() {
  static_assert(std::is_same<Ret, arma::mat>::value,
                "Fun must be of type arma::mat (const arma::mat&)");
}

template <typename Fun>
  requires requires(Fun fun) {
    { fun(std::declval<arma::mat const &>()) } -> std::same_as<arma::mat>;
  }
std::tuple<WeightedClusters, double>
GraphCut::partitionGraph(std::uint8_t nClusters, float weightModule,
                         std::uint16_t nTries, Fun graphFun) const {
  constexpr unsigned minClusterTotalWeight = 3;

  auto const &firstAdjacency = adjacencies[0];
  assert(std::ranges::all_of(adjacencies | std::views::drop(1),
                             [&](auto const &adjacency) {
                               return adjacency.n_rows == firstAdjacency.n_rows;
                             }));

  checkGraphFunCallable<std::decay_t<Fun>>();
  nClusters =
      std::min(nClusters, static_cast<std::uint8_t>(std::min(
                              firstAdjacency.n_rows / minClusterTotalWeight,
                              static_cast<unsigned long long>(
                                  std::numeric_limits<std::uint8_t>::max()))));
  if (nClusters < 2) {
    if (auto init = std::get_if<WeightedClusters>(&initialClusters))
      return std::tuple{*init, NAN};
    else
      return std::tuple{WeightedClusters(firstAdjacency.n_rows, 1), NAN};
  }

  auto graphs = adjacencies | std::views::transform([=](auto const &adjacency) {
                  return graphFun(adjacency);
                }) |
                std::views::as_rvalue | std::ranges::to<std::vector>();

  WeightedClusters weights = [&, this] {
    if (auto init = std::get_if<WeightedClusters>(&initialClusters))
      return *init;
    else {
      const std::size_t nBases = firstAdjacency.n_rows;
      WeightedClusters weights(nBases, nClusters, false);
      std::mt19937 randomGen(std::random_device{}());
      // std::uniform_int_distribution<std::size_t> clusterAssigner(0,
      //                                                           nClusters -
      //                                                           1);

      std::vector<std::size_t> usableBases(nClusters, 0);
      std::size_t clusterUsableBasesIndex = 0;
      std::size_t maxUsableBases = nBases;
      const std::size_t usableBasesModule =
          std::max(nBases / (nClusters * 5),
                   static_cast<std::size_t>(minClusterTotalWeight));

      auto incrementUsableBases = [&](std::size_t index) {
        auto r = usableBases | std::views::take(index);
        usableBases[index] = std::min(
            usableBases[index] + usableBasesModule,
            nBases - std::accumulate(std::ranges::begin(r), std::ranges::end(r),
                                     std::size_t(0)));
      };

      usableBases[0] = usableBasesModule;
      auto bestScore = std::numeric_limits<double>::infinity();
      WeightedClusters temporaryWeights(nBases, nClusters, false);
      std::vector<std::size_t> baseIndices(nBases);
      std::iota(std::begin(baseIndices), std::end(baseIndices), std::size_t(0));

      for (;;) {
        if (usableBases[clusterUsableBasesIndex] == 0) {
          auto r = usableBases | std::views::take(clusterUsableBasesIndex + 1);

          maxUsableBases =
              nBases - std::accumulate(std::ranges::begin(r),
                                       std::ranges::end(r), std::size_t(0));

          assert(maxUsableBases <= nBases);
          if (clusterUsableBasesIndex ==
              static_cast<std::size_t>(nClusters - 1)) {
            usableBases[clusterUsableBasesIndex] = maxUsableBases;
            maxUsableBases = 0;
          } else {
            usableBases[clusterUsableBasesIndex] =
                std::min(usableBasesModule, maxUsableBases);

            assert(maxUsableBases >= usableBases[clusterUsableBasesIndex]);
            maxUsableBases -= usableBases[clusterUsableBasesIndex];
          }
        }

        if (maxUsableBases > 0) {
          ++clusterUsableBasesIndex;
          assert(clusterUsableBasesIndex < nClusters);
          continue;
        }

        if ((nClusters > 0 and clusterUsableBasesIndex <
                                   static_cast<std::size_t>(nClusters - 1)) or
            usableBases.back() == 0) {
          for (;;) {
            if (usableBases[clusterUsableBasesIndex] <= usableBasesModule) {
              usableBases[clusterUsableBasesIndex] = 0;
              --clusterUsableBasesIndex;
            } else {
              if (usableBases[clusterUsableBasesIndex] >=
                  nBases - usableBasesModule) {
                usableBases[clusterUsableBasesIndex] = 0;
                if (clusterUsableBasesIndex > 0) {
                  auto lastUsableBasesIndex = clusterUsableBasesIndex - 1;

                  assert(lastUsableBasesIndex < nClusters);
                  incrementUsableBases(lastUsableBasesIndex);

                  assert(std::accumulate(std::begin(usableBases),
                                         std::end(usableBases),
                                         std::size_t(0)) <= nBases);
                }
                break;
              } else {
                if (clusterUsableBasesIndex > 0) {
                  assert(clusterUsableBasesIndex - 1 < nClusters);

                  incrementUsableBases(clusterUsableBasesIndex - 1);
                  usableBases[clusterUsableBasesIndex] = 0;
                  assert(std::accumulate(std::begin(usableBases),
                                         std::end(usableBases),
                                         std::size_t(0)) <= nBases);
                  break;
                } else {
                  usableBases[0] = 0;
                  break;
                }
              }
            }
          }

          if (usableBases[0] == 0)
            break;
          else
            continue;
        }

        assert(std::ranges::all_of(
            usableBases, [&](std::size_t bases) { return bases <= nBases; }));
        assert(std::accumulate(std::begin(usableBases), std::end(usableBases),
                               std::size_t(0)) == nBases);
        assert(std::ranges::none_of(
            usableBases, [](std::size_t nBases) { return nBases == 0; }));
        for (std::uint16_t trialIndex = 0; trialIndex < nTries; ++trialIndex) {
          temporaryWeights.removeWeights();

          std::ranges::shuffle(baseIndices, randomGen);

          auto baseIndicesIter = std::ranges::begin(baseIndices);
          auto usableBasesIter = std::ranges::begin(usableBases);
          for (unsigned clusterIndex = 0; clusterIndex < nClusters;
               ++clusterIndex, ++usableBasesIter) {
            assert(usableBasesIter < std::ranges::end(usableBases));
            assert(baseIndicesIter < std::ranges::end(baseIndices));

            const auto currentUsableBases =
                static_cast<std::ptrdiff_t>(*usableBasesIter);
            auto &&currentCluster = temporaryWeights.cluster(clusterIndex);
            const auto currentBaseIndicesEnd =
                std::ranges::next(baseIndicesIter, currentUsableBases);
            assert(
                static_cast<std::size_t>(std::ranges::distance(
                    std::ranges::begin(baseIndices), currentBaseIndicesEnd)) ==
                std::accumulate(std::ranges::begin(usableBases),
                                std::ranges::next(usableBasesIter),
                                std::size_t(0)));
            assert(currentBaseIndicesEnd <= std::ranges::end(baseIndices));

            std::ranges::for_each(baseIndicesIter, currentBaseIndicesEnd,
                                  [&](std::size_t baseIndex) {
                                    currentCluster[baseIndex] = 1.f;
                                  });

            baseIndicesIter = currentBaseIndicesEnd;
          }

          double score = std::ranges::fold_left(
              graphs | std::views::transform([&](auto const &graph) {
                return calculateCutScore(graph, temporaryWeights);
              }),
              0., std::plus{});
          assert(not std::isnan(score));

          if (score < bestScore and not std::isinf(score)) {
            bestScore = score;
            weights = temporaryWeights;
          }
        }

        usableBases.back() = 0;
        incrementUsableBases(nClusters - 2);
      }
      return weights;
    }
  }();

  auto weightsClusters = weights.clusters();
  static_assert(not std::is_const<decltype(weightsClusters)>::value, "!");
  float currentWeightChange = weightModule;
  double bestScore = NAN;

  for (const auto clustersEnd = std::ranges::end(weightsClusters);;) {
    bestScore = std::ranges::fold_left(
        graphs | std::views::transform([&](auto const &graph) {
          return calculateCutScore(graph, weights);
        }),
        0., std::plus{});

    if (std::isnan(bestScore)) {
      std::stringstream graphs_content;

      auto graphs_to_dump =
          graphs | std::views::transform([&](auto const &graph) {
            auto score = calculateCutScore(graph, weights);
            return std::pair(&graph, score);
          }) |
          std::views::filter(
              [](auto const &pair) { return std::isnan(std::get<1>(pair)); }) |
          std::views::transform([](auto &&pair) { return std::get<0>(pair); });
      (*std::ranges::begin(graphs_to_dump))->print(graphs_content);
      for (auto graph_to_dump : graphs_to_dump | std::views::drop(1)) {
        graphs_content << "\n\n";
        graph_to_dump->print(graphs_content);
      }

      logger::error("Invalid bestScore. Graphs content:\n{}",
                    graphs_content.str());
      std::abort();
    }

    auto bestFromClusterIter = clustersEnd;
    auto bestToClusterIter = bestFromClusterIter;

    // FIXME GCC warns about this that could be used uninitialized
    // (theoretically, it cannot happens...). Let's init it...
    auto bestBaseIndex = std::numeric_limits<std::size_t>::max();

    for (auto fromClusterIter = std::ranges::begin(weightsClusters);
         fromClusterIter < clustersEnd; ++fromClusterIter) {
      auto &&fromCluster = *fromClusterIter;

      const auto fromClusterEnd = std::ranges::end(fromCluster);
      if (std::accumulate(std::ranges::begin(fromCluster), fromClusterEnd,
                          0.f) <= 3)
        continue;

      for (auto toClusterIter = std::ranges::begin(weightsClusters);
           toClusterIter < clustersEnd; ++toClusterIter) {
        if (fromClusterIter == toClusterIter)
          continue;

        auto &&toCluster = *toClusterIter;

        auto fromClusterCurrentIter = std::ranges::begin(fromCluster);
        auto toClusterCurrentIter = std::ranges::begin(toCluster);
        std::size_t baseIndex = 0;
        for (; fromClusterCurrentIter < fromClusterEnd;
             ++fromClusterCurrentIter, ++toClusterCurrentIter, ++baseIndex) {
          auto &fromClusterWeight = *fromClusterCurrentIter;
          if (fromClusterWeight == 0.f or
              fromClusterWeight < currentWeightChange)
            continue;
          auto &toClusterWeight = *toClusterCurrentIter;
          if (toClusterWeight > 1.f - currentWeightChange)
            continue;

          const auto fromClusterOldWeight = fromClusterWeight;
          const auto toClusterOldWeight = toClusterWeight;

          fromClusterWeight -= currentWeightChange;
          toClusterWeight += currentWeightChange;

          double currentScore = std::ranges::fold_left(
              graphs | std::views::transform([&](auto const &graph) {
                return calculateCutScore(graph, weights);
              }),
              0., std::plus{});
          assert(not std::isnan(currentScore));

          if (currentScore < bestScore) {
            bestScore = currentScore;
            bestFromClusterIter = fromClusterIter;
            bestToClusterIter = toClusterIter;
            bestBaseIndex = baseIndex;
          }

          fromClusterWeight = fromClusterOldWeight;
          toClusterWeight = toClusterOldWeight;
        }
      }
    }

    if (bestFromClusterIter == clustersEnd) {
      currentWeightChange += weightModule;
      if (currentWeightChange > 1.f)
        break;

      continue;
    }

    {
      auto &&bestFromCluster = *bestFromClusterIter;
      auto &&bestToCluster = *bestToClusterIter;

      // FIXME this should not assert, but it is related to the suspicious
      // warning of GCC
      assert(bestBaseIndex != std::numeric_limits<std::size_t>::max());
      bestFromCluster[bestBaseIndex] -= currentWeightChange;
      bestToCluster[bestBaseIndex] += currentWeightChange;

      currentWeightChange = weightModule;
    }
  }

  return std::tuple{weights, bestScore};
}

template <typename Fun>
HardClusters GraphCut::partitionGraphHard(std::uint8_t nClusters,
                                          Fun graphFun) const {
  checkGraphFunCallable<std::decay_t<Fun>>();
  auto const &firstAdjacency = adjacencies[0];
  assert(std::ranges::all_of(adjacencies | std::views::drop(1),
                             [&](auto const &adjacency) {
                               return adjacency.n_rows == firstAdjacency.n_rows;
                             }));

  auto graphs = adjacencies | std::views::transform([&](auto const &adjacency) {
                  return graphFun(adjacency);
                }) |
                std::ranges::to<std::vector>();

  HardClusters hardClusters = [&, this] {
    if (auto init = std::get_if<HardClusters>(&initialClusters))
      return *init;
    else {
      HardClusters hardClusters(firstAdjacency.n_rows, nClusters);
      std::mt19937 randomGen(std::random_device{}());
      std::uniform_int_distribution<std::uint8_t> clusterAssigner(0, nClusters -
                                                                         1);
      auto hardClustersWrapper = hardClusters.clusters();
      do {
        hardClusters.clear();

        for (unsigned baseIndex = 0; baseIndex < firstAdjacency.n_rows;
             ++baseIndex)
          hardClusters[baseIndex].set(clusterAssigner(randomGen));
      } while (std::ranges::any_of(hardClustersWrapper,
                                   [](const auto &cluster) {
                                     return std::ranges::count(cluster, true) <
                                            3;
                                   }) or
               std::ranges::any_of(graphs, [&](auto const &graph) {
                 return std::isinf(calculateCutScore(graph, hardClusters));
               }));

      return hardClusters;
    }
  }();

  auto clusters = hardClusters.clusters();

  for (const auto clustersEnd = std::ranges::end(clusters);;) {
    double bestScore = std::ranges::fold_left(
        graphs | std::views::transform([&](auto const &graph) {
          return calculateCutScore(graph, hardClusters);
        }),
        0., std::plus{});
    if (std::isnan(bestScore)) {
      std::stringstream graphs_content;

      auto graphs_to_dump =
          graphs | std::views::transform([&](auto const &graph) {
            auto score = calculateCutScore(graph, hardClusters);
            return std::pair(&graph, score);
          }) |
          std::views::filter(
              [](auto const &pair) { return std::isnan(std::get<1>(pair)); }) |
          std::views::transform([](auto &&pair) { return std::get<0>(pair); });
      (*std::ranges::begin(graphs_to_dump))->print(graphs_content);
      for (auto graph_to_dump : graphs_to_dump | std::views::drop(1)) {
        graphs_content << "\n\n";
        graph_to_dump->print(graphs_content);
      }

      logger::error("Invalid bestScore. Graphs content:\n{}",
                    graphs_content.str());
      std::abort();
    }

    auto bestFromClusterIter = clustersEnd;
    auto bestToClusterIter = bestFromClusterIter;
    // FIXME GCC warns about this that could be used uninitialized
    // (theoretically, it cannot happens...). Let's init it...
    auto bestBaseIndex = std::numeric_limits<std::size_t>::max();

    for (auto fromClusterIter = std::ranges::begin(clusters);
         fromClusterIter < clustersEnd; ++fromClusterIter) {
      auto fromCluster = *fromClusterIter;
      if (std::ranges::count(fromCluster, true) <= 3)
        continue;

      const auto fromClusterEnd = std::ranges::end(fromCluster);

      for (auto toClusterIter = std::ranges::begin(clusters);
           toClusterIter < clustersEnd; ++toClusterIter) {
        if (fromClusterIter == toClusterIter)
          continue;

        auto toCluster = *toClusterIter;

        auto fromClusterCurrentIter = std::ranges::begin(fromCluster);
        auto toClusterCurrentIter = std::ranges::begin(toCluster);
        std::size_t baseIndex = 0;
        for (; fromClusterCurrentIter < fromClusterEnd;
             ++fromClusterCurrentIter, ++toClusterCurrentIter, ++baseIndex) {
          auto fromClusterCurrent = *fromClusterCurrentIter;
          if (not fromClusterCurrent)
            continue;

          auto toClusterCurrent = *toClusterCurrentIter;
          if (toClusterCurrent)
            continue;

          fromClusterCurrent.swap(toClusterCurrent);

          double currentScore = std::ranges::fold_left(
              graphs | std::views::transform([&](auto const &graph) {
                return calculateCutScore(graph, hardClusters);
              }),
              0., std::plus{});
          assert(not std::isnan(currentScore));

          if (not std::isinf(currentScore) and currentScore < bestScore) {
            bestScore = currentScore;
            bestFromClusterIter = fromClusterIter;
            bestToClusterIter = toClusterIter;
            bestBaseIndex = baseIndex;
          }

          fromClusterCurrent.swap(toClusterCurrent);
        }
      }
    }

    if (bestFromClusterIter == clustersEnd)
      return hardClusters;

    {
      auto bestFromCluster = *bestFromClusterIter;
      auto bestToCluster = *bestToClusterIter;

      // FIXME this should not assert, but it is related to the suspicious
      // warning of GCC
      assert(bestBaseIndex != std::numeric_limits<std::size_t>::max());
      bestFromCluster[bestBaseIndex].swap(bestToCluster[bestBaseIndex]);
    }
  }
}

template <typename Clusters>
double GraphCut::calculateCutScore(const arma::mat &graph,
                                   const Clusters &clusters) const {
  static_assert(is_clusters_v<Clusters>);
  double score = 0;
  auto complement = clusters.complement();
  auto clustersWrapper = clusters.clusters();
  for (const auto &cluster : clustersWrapper) {
    double divisor = calculateCut(graph, cluster, cluster);
    if (divisor == 0.)
      return std::numeric_limits<double>::infinity();

    score += calculateCut(graph, cluster, complement.cluster(cluster.index())) /
             divisor;
    assert(not std::isnan(score));
  }

  return score;
}

template <typename ClusterA, typename ClusterB>
double GraphCut::calculateCut(const arma::mat &graph, const ClusterA &clusterA,
                              const ClusterB &clusterB) {
  static_assert(is_cluster_wrapper_v<ClusterA>);
  static_assert(is_cluster_wrapper_v<ClusterB>);
  double score = 0;
  using value_type = typename ClusterA::iterator::value_type;

  static constexpr bool valueIsIntegral = std::is_integral_v<value_type>;
  static_assert(valueIsIntegral ==
                std::is_integral_v<typename ClusterB::iterator::value_type>);

  std::size_t indexA = 0;
  for (auto &&valueA : clusterA) {
    const bool skip = [&] {
      if constexpr (valueIsIntegral)
        return not valueA;
      else
        return valueA == 0;
    }();

    if (skip) {
      ++indexA;
      continue;
    }

    const auto &colA = graph.col(indexA);
    std::size_t indexB = 0;
    for (auto &&valueB : clusterB) {
      if constexpr (valueIsIntegral) {
        if (valueB)
          score += colA(indexB);
      } else
        score += colA(indexB) * valueA * valueB;
      ++indexB;
    }

    ++indexA;
  }

  assert(not std::isnan(score));
  return score;
}

template <typename Clusters>
void GraphCut::setInitialClusters(Clusters &&clusters) {
  using clusters_type = std::decay_t<Clusters>;
  static_assert(std::is_same_v<clusters_type, HardClusters> or
                std::is_same_v<clusters_type, WeightedClusters>);

  initialClusters.emplace<clusters_type>(std::forward<Clusters>(clusters));
}

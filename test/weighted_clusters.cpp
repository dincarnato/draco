#include "weighted_clusters.hpp"

#include <range/v3/view/reverse.hpp>

#include <array>
#include <cassert>
#include <utility>

static constexpr std::ptrdiff_t nClusters = 3;
static constexpr std::ptrdiff_t nElements = 10;
static constexpr std::array<float, nClusters * nElements> allWeights{
    0.50f, 0.30f, 0.20f, 0.51f, 0.31f, 0.18f, 0.52f, 0.32f, 0.16f, 0.53f,
    0.33f, 0.14f, 0.54f, 0.34f, 0.12f, 0.55f, 0.35f, 0.10f, 0.56f, 0.36f,
    0.08f, 0.57f, 0.37f, 0.06f, 0.58f, 0.38f, 0.04f, 0.59f, 0.39f, 0.02f,
};

template <bool post_increment, bool iter_forward, typename T>
static void
test_weighted_clusters(T&& weightedClusters) {
  static_assert(std::is_same_v<std::decay_t<T>, WeightedClusters>);

  auto weightsRange = [&] {
    if constexpr (iter_forward)
      return allWeights;
    else
      return allWeights | ranges::view::reverse;
  }();

  auto getRange = [](auto&& range) {
    if constexpr (iter_forward)
      return range;
    else
      return range | ranges::view::reverse;
  };

  auto getWeightsIter = [](auto& weightsRange, std::ptrdiff_t index) {
    if constexpr (not iter_forward) {
      assert(index < nClusters);
      index = nClusters - index - 1;
    }
    return ranges::next(ranges::begin(weightsRange), index);
  };

  {
    std::ptrdiff_t clusterIndex = [] {
      if constexpr (iter_forward)
        return std::ptrdiff_t(0);
      else
        return std::ptrdiff_t(nClusters - 1);
    }();
    auto&& clusters = weightedClusters.clusters();
    assert(clusters.size() == nClusters);
    auto clustersRange = getRange(clusters);
    for (
        auto clusterIter = ranges::begin(clustersRange);
        clusterIter != ranges::end(clustersRange); [&] {
          if constexpr (post_increment)
            clusterIter++;
          else
            ++clusterIter;
        }(),
             [&] {
               if constexpr (iter_forward)
                 ++clusterIndex;
               else
                 --clusterIndex;
             }()) {

      assert(clusterIter >= ranges::begin(clustersRange));
      assert(clusterIter < ranges::end(clustersRange));

      auto&& cluster = *clusterIter;
      assert(cluster.size() == nElements);
      auto weightsIter = getWeightsIter(weightsRange, clusterIndex);
      auto clusterRange = getRange(cluster);
      for (auto elementIter = ranges::begin(clusterRange);
           elementIter != ranges::end(clusterRange); [&] {
             if constexpr (post_increment)
               elementIter++;
             else
               ++elementIter;
           }()) {
        assert(weightsIter >= ranges::begin(weightsRange) and
               weightsIter < ranges::end(weightsRange));
        assert(elementIter >= ranges::begin(clusterRange));
        assert(elementIter < ranges::end(clusterRange));

        auto&& element = *elementIter;
        assert(element == *weightsIter);

        if (ranges::distance(weightsIter, ranges::end(weightsRange)) >=
            nClusters)
          ranges::advance(weightsIter, nClusters);
        else
          weightsIter = ranges::end(weightsRange);
      }
    }
  }

  for (std::ptrdiff_t clusterIndex =
           [] {
             if constexpr (iter_forward)
               return std::ptrdiff_t(0);
             else
               return std::ptrdiff_t(nClusters - 1);
           }();
       [&] {
         if constexpr (iter_forward)
           return clusterIndex < nClusters;
         else
           return clusterIndex >= 0;
       }();
       [&] {
         if constexpr (iter_forward)
           ++clusterIndex;
         else
           --clusterIndex;
       }()) {
    auto&& cluster =
        weightedClusters.cluster(static_cast<std::size_t>(clusterIndex));
    auto weightsIter = getWeightsIter(weightsRange, clusterIndex);

    assert(cluster.size() == nElements);
    auto clusterRange = getRange(cluster);
    for (auto valueIter = ranges::begin(clusterRange);
         valueIter != ranges::end(clusterRange); [&] {
           if constexpr (post_increment)
             valueIter++;
           else
             ++valueIter;
         }()) {
      assert(weightsIter >= ranges::begin(weightsRange) and
             weightsIter < ranges::end(weightsRange));
      assert(valueIter >= ranges::begin(clusterRange));
      assert(valueIter < ranges::end(clusterRange));

      auto&& value = *valueIter;
      assert(value == *weightsIter);

      if (ranges::distance(weightsIter, ranges::end(weightsRange)) >= nClusters)
        ranges::advance(weightsIter, nClusters);
      else
        weightsIter = ranges::end(weightsRange);
    }
  }

  {
    auto clustersWrapper = weightedClusters.clusters();
    auto firstClusterIter = ranges::begin(clustersWrapper);

    for (std::ptrdiff_t clusterIndex =
             [] {
               if constexpr (iter_forward)
                 return std::ptrdiff_t(0);
               else
                 return std::ptrdiff_t(nClusters - 1);
             }();
         [&] {
           if constexpr (iter_forward)
             return clusterIndex < nClusters;
           else
             return clusterIndex >= 0;
         }();
         [&] {
           if constexpr (iter_forward)
             ++clusterIndex;
           else
             --clusterIndex;
         }()) {
      auto&& cluster = firstClusterIter[clusterIndex];
      auto clusterFirstElementIter = ranges::begin(cluster);

      for (std::ptrdiff_t elementIndex =
               [] {
                 if constexpr (iter_forward)
                   return std::ptrdiff_t(0);
                 else
                   return std::ptrdiff_t(nElements - 1);
               }();
           [&] {
             if constexpr (iter_forward)
               return elementIndex < nElements;
             else
               return elementIndex >= 0;
           }();
           [&] {
             if constexpr (iter_forward)
               ++elementIndex;
             else
               --elementIndex;
           }()) {

        auto&& value = clusterFirstElementIter[elementIndex];
        assert(value == allWeights[static_cast<std::size_t>(
                            elementIndex * nClusters + clusterIndex)]);
      }
    }
  }

  {
    auto firstElementIter = ranges::begin(weightedClusters);

    for (std::ptrdiff_t elementIndex =
             [] {
               if constexpr (iter_forward)
                 return std::ptrdiff_t(0);
               else
                 return std::ptrdiff_t(nElements - 1);
             }();
         [&] {
           if constexpr (iter_forward)
             return elementIndex < nElements;
           else
             return elementIndex >= 0;
         }();
         [&] {
           if constexpr (iter_forward)
             ++elementIndex;
           else
             --elementIndex;
         }()) {

      auto&& elementSpan = firstElementIter[elementIndex];
      auto elementSpanFirstClusterIter = ranges::begin(elementSpan);

      for (std::ptrdiff_t clusterIndex =
               [] {
                 if constexpr (iter_forward)
                   return std::ptrdiff_t(0);
                 else
                   return std::ptrdiff_t(nClusters - 1);
               }();
           [&] {
             if constexpr (iter_forward)
               return clusterIndex < nClusters;
             else
               return clusterIndex >= 0;
           }();
           [&] {
             if constexpr (iter_forward)
               ++clusterIndex;
             else
               --clusterIndex;
           }()) {

        auto&& value =
            elementSpanFirstClusterIter[static_cast<std::size_t>(clusterIndex)];
        assert(value == allWeights[static_cast<std::size_t>(
                            elementIndex * nClusters + clusterIndex)]);
      }
    }
  }

  for (unsigned elementIndex = 0; elementIndex < nElements; ++elementIndex) {
    auto&& elementSpan = weightedClusters[elementIndex];

    assert(ranges::distance(ranges::begin(elementSpan),
                            ranges::end(elementSpan)) == nClusters);
    assert(elementSpan.span_size() == nClusters);
    for (unsigned clusterIndex = 0; clusterIndex < nClusters; ++clusterIndex)
      assert(elementSpan[static_cast<std::size_t>(clusterIndex)] ==
             allWeights[static_cast<std::size_t>(elementIndex * nClusters +
                                                 clusterIndex)]);
  }
}

int
main() {
  {
    WeightedClusters weightedClusters;
    assert(weightedClusters.getClustersSize() == 0);
    assert(weightedClusters.getElementsSize() == 0);
  }

  {
    WeightedClusters weightedClusters(3, 3, false);
    assert(weightedClusters.getClustersSize() == 3);
    assert(weightedClusters.getElementsSize() == 3);

    for (auto&& cluster : weightedClusters) {
      for (auto&& value : cluster)
        assert(value == 0);
    }
  }

  WeightedClusters weightedClusters(nElements, nClusters);

  assert(weightedClusters.getClustersSize() == nClusters);
  assert(weightedClusters.getElementsSize() == nElements);

  {
    auto iter = ranges::begin(weightedClusters);
    auto allWeightsIter = ranges::begin(allWeights);

    for (unsigned elementIndex = 0; elementIndex < nElements; ++elementIndex) {
      auto elementWrapper = *iter++;
      auto&& elementIter = ranges::begin(elementWrapper);
      for (unsigned clusterIndex = 0; clusterIndex < nClusters; ++clusterIndex)
        *elementIter++ = *allWeightsIter++;
    }
  }

  test_weighted_clusters<false, true>(weightedClusters);
  test_weighted_clusters<false, true>(std::as_const(weightedClusters));
  test_weighted_clusters<false, false>(weightedClusters);
  test_weighted_clusters<false, false>(std::as_const(weightedClusters));
  test_weighted_clusters<true, true>(weightedClusters);
  test_weighted_clusters<true, true>(std::as_const(weightedClusters));
  test_weighted_clusters<true, false>(weightedClusters);
  test_weighted_clusters<true, false>(std::as_const(weightedClusters));
}

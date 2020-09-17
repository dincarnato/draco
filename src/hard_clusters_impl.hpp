#pragma once

#include "hard_clusters.hpp"

#include <algorithm>
#include <cassert>
#include <limits>

template <typename Type>
HardClustersBase<Type>::HardClustersBase(std::size_t elements,
                                         Type clusters) noexcept
    : nClusters(clusters),
      _clusters(elements, std::numeric_limits<Type>::max()) {
  assert(clusters < std::numeric_limits<Type>::max());
}

template <typename Type>
HardClustersBase<Type>::HardClustersBase(std::size_t elements, Type clusters,
                                         Type defaultCluster) noexcept
    : nClusters(clusters), _clusters(elements, defaultCluster) {
  assert(clusters < std::numeric_limits<Type>::max());
  assert(defaultCluster < clusters);
}

template <typename Type>
auto
HardClustersBase<Type>::cluster(index_type index) -> cluster_wrapper_type {
  assert(index < nClusters);
  return {*this, index};
}

template <typename Type>
auto
HardClustersBase<Type>::cluster(index_type index) const
    -> const_cluster_wrapper_type {
  assert(index < nClusters);
  return {*this, index};
}

template <typename Type>
void
HardClustersBase<Type>::clear() {
  std::fill(std::begin(_clusters), std::end(_clusters),
            std::numeric_limits<Type>::max());
}

template <typename Type>
auto HardClustersBase<Type>::operator[](std::size_t index) -> element_wrapper {
  return {*this, index};
}

template <typename Type>
auto HardClustersBase<Type>::operator[](std::size_t index) const
    -> const_element_wrapper {
  return {*this, index};
}

template <typename Type>
auto
HardClustersBase<Type>::clusters() -> clusters_wrapper_type {
  return clusters_wrapper_type{*this};
}

template <typename Type>
auto
HardClustersBase<Type>::clusters() const -> const_clusters_wrapper_type {
  return const_clusters_wrapper_type{*this};
}

template <typename Type>
auto
HardClustersBase<Type>::complement() -> complement_type {
  return complement_type{*this};
}

template <typename Type>
auto
HardClustersBase<Type>::complement() const -> const_complement_type {
  return const_complement_type{*this};
}

template <typename Type>
std::size_t
HardClustersBase<Type>::getElementsSize() const {
  return _clusters.size();
}

template <typename Type>
auto
HardClustersBase<Type>::getClustersSize() const -> index_type {
  return nClusters;
}

template <typename Type>
auto
HardClustersBase<Type>::begin() -> iterator {
  return std::begin(_clusters);
}

template <typename Type>
auto
HardClustersBase<Type>::end() -> iterator {
  return std::end(_clusters);
}

template <typename Type>
auto
HardClustersBase<Type>::begin() const -> const_iterator {
  return std::begin(_clusters);
}

template <typename Type>
auto
HardClustersBase<Type>::end() const -> const_iterator {
  return std::end(_clusters);
}

template <typename Type>
void
HardClustersBase<Type>::extendMinorCluster(std::size_t maxDistance) {
  /*
  const index_type minorCluster = [&] {
    std::vector<std::size_t> counts(nClusters, 0);
    for (index_type cluster : _clusters)
      ++counts[cluster];
    return std::distance(
        std::begin(counts),
        std::min_element(std::begin(counts), std::end(counts)));
  }();
  */

  struct Hole {
    iterator begin;
    iterator end;
    std::size_t size;

    bool
    operator<(const Hole& other) const {
      return (size < other.size) or
             (size == other.size and begin < other.begin);
    }
  };
  std::vector<Hole> holes;

  const auto beginClusters = std::begin(_clusters);
  const auto endClusters = std::end(_clusters);
  {
    auto spanStart = beginClusters;
    for (auto spanEnd = std::next(spanStart); spanStart != endClusters;
         ++spanEnd) {
      if (spanEnd == endClusters or *spanEnd != *spanStart) {
        if (std::size_t holeSize = std::distance(spanStart, spanEnd);
            holeSize <= maxDistance)
          holes.push_back({spanStart, spanEnd, holeSize});

        spanStart = spanEnd;
      }
    }
  }

  std::sort(std::begin(holes), std::end(holes));

  const auto beginHoles = std::begin(holes);
  auto endHoles = std::end(holes);

  for (auto currentHole = beginHoles; currentHole < endHoles;) {
    auto nextHole = currentHole;
    for (;;) {
      auto nextNextHole = std::next(nextHole);
      if (nextNextHole != endHoles and nextHole->size == nextNextHole->size and
          nextHole->end == nextNextHole->begin)
        nextHole = nextNextHole;
      else
        break;
    }

    auto lastHole = nextHole++;
    const auto& currentHoleBegin = currentHole->begin;
    const auto& currentHoleEnd = lastHole->end;

    auto updateModifiedHole = [&](iterator elementIter) {
      auto holeIter = std::find_if(nextHole, endHoles, [=](const auto& hole) {
        return hole.begin == elementIter or hole.end == elementIter;
      });

      if (holeIter != endHoles) {
        auto& hole = *holeIter;
        if (hole.begin == elementIter)
          hole.begin = currentHole->begin;
        else
          hole.end = lastHole->end;

        hole.size = std::distance(hole.begin, hole.end);
      }

      return holeIter;
    };

    auto spanBeforeSize = [&] {
      if (currentHoleBegin == beginClusters)
        return static_cast<std::size_t>(0);
      else {
        auto prevElementIter = std::prev(currentHoleBegin);
        const auto& prevElement = *prevElementIter;

        auto prevSpanBegin = prevElementIter;
        for (;;) {
          if (prevSpanBegin < beginClusters or *prevSpanBegin != prevElement)
            break;

          --prevSpanBegin;
        }
        ++prevSpanBegin;

        return static_cast<std::size_t>(
            std::distance(prevSpanBegin, currentHoleBegin));
      }
    }();

    auto spanAfterSize = [&] {
      auto nextElementIter = currentHoleEnd;
      if (nextElementIter == endClusters)
        return static_cast<std::size_t>(0);
      else {
        const auto& nextElement = *nextElementIter;

        auto nextSpanEnd = std::next(nextElementIter);
        for (;;) {
          if (nextSpanEnd == endClusters or *nextSpanEnd != nextElement)
            break;

          ++nextSpanEnd;
        }

        return static_cast<std::size_t>(
            std::distance(nextElementIter, nextSpanEnd));
      }
    }();

    bool changedHoles = false;
    auto spanToCopyFrom = [&] {
      auto prevCurrentHoleBegin = std::prev(currentHoleBegin);

      if (spanBeforeSize != 0 and spanAfterSize != 0) {
        const auto& clusterBefore = *prevCurrentHoleBegin;
        const auto& clusterAfter = *currentHoleEnd;

        if (clusterBefore == clusterAfter) {
          if (auto holeAfter = updateModifiedHole(currentHoleEnd);
              holeAfter != endHoles) {

            changedHoles = true;
            if (auto holeBefore = updateModifiedHole(prevCurrentHoleBegin);
                holeBefore != endHoles) {
              endHoles = holes.erase(holeAfter);
            }
          } else if (updateModifiedHole(prevCurrentHoleBegin) != endHoles)
            changedHoles = true;

          return currentHoleEnd;
        } else {
          if (spanAfterSize < spanBeforeSize) {
            /*
             * We are choosing to expand the smaller hole. Is it a good idea?
             */
            if (updateModifiedHole(currentHoleEnd) != endHoles)
              changedHoles = true;

            return currentHoleEnd;
          } else {
            /*
             * In case they have the same size, we are choosing the first span
             * randomly, we don't have a specific way to handle this situation
             * (yet)
             */
            if (updateModifiedHole(prevCurrentHoleBegin) != endHoles)
              changedHoles = true;

            return prevCurrentHoleBegin;
          }
        }
      } else if (spanBeforeSize != 0) {
        if (updateModifiedHole(prevCurrentHoleBegin) != endHoles)
          changedHoles = true;

        return prevCurrentHoleBegin;
      } else {
        if (updateModifiedHole(currentHoleEnd) != endHoles)
          changedHoles = true;

        return currentHoleEnd;
      }
    }();

    std::fill(currentHoleBegin, currentHoleEnd, *spanToCopyFrom);
    if (changedHoles) {
      std::sort(nextHole, endHoles);
      endHoles = std::find_if(nextHole, endHoles, [=](const auto& hole) {
        return hole.size > maxDistance;
      });
    }
    currentHole = nextHole;
  }
}

template <typename Type>
void
HardClustersBase<Type>::removeSmallSpans(std::size_t minSize) {
  const auto beginClusters = std::begin(_clusters);
  const auto endClusters = std::end(_clusters);
  for (auto clusterIter = beginClusters; clusterIter < endClusters;) {
    const auto spanEnd =
        std::find_if_not(std::next(clusterIter), endClusters,
                         [clusterIndex = *clusterIter](index_type cluster) {
                           return cluster == clusterIndex;
                         });
    if (static_cast<std::size_t>(std::distance(clusterIter, spanEnd)) <
        minSize) {
      // TODO better approach than extending the first cluster (if possible)
      if (clusterIter != beginClusters)
        std::fill(clusterIter, spanEnd, *std::prev(clusterIter));
      else if (spanEnd != endClusters)
        std::fill(clusterIter, spanEnd, *spanEnd);
      else
        throw std::runtime_error("cannot remove the only span available");
    }

    clusterIter = spanEnd;
  }
}

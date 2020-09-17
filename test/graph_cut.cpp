#include "graph_cut.hpp"

#include <range/v3/algorithm.hpp>

#include <algorithm>
#include <armadillo>
#include <array>
#include <vector>

static const arma::mat adjacency{{0, 0, 0, 1, 4, 0, 3, 0, 3, 1, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 4, 1, 0, 0, 2, 3, 0, 1},
                                 {0, 0, 0, 2, 0, 2, 0, 4, 0, 0, 0, 4, 1},
                                 {0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 4, 0},
                                 {0, 0, 0, 0, 0, 0, 4, 0, 4, 4, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 3, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

static const std::array<std::vector<int>, 3> expectedResults{
    std::vector{0, 4, 6, 8, 9}, std::vector{1, 5, 10},
    std::vector{2, 3, 7, 11, 12}};

int
main() {
  namespace rng = ::ranges;

  GraphCut graphCut(arma::symmatu(adjacency));
  auto results = graphCut.run(3, GraphCut::hard);
  if (results.getClustersSize() == 0 or results.getElementsSize() == 0)
    return -1;

  {
    auto clustersWrapper = results.clusters();
    if (rng::any_of(rng::begin(clustersWrapper), rng::end(clustersWrapper),
                    [](const auto& cluster) {
                      return rng::distance(rng::begin(cluster),
                                           rng::end(cluster)) == 0;
                    }))
      return -1;
  }

  std::array<std::vector<int>, 3> parsedResults;
  auto const elements_size = static_cast<int>(results.getElementsSize());
  for (int elementIndex = 0; elementIndex < elements_size; ++elementIndex) {
    const auto assignment =
        results[static_cast<std::size_t>(elementIndex)].get();
    parsedResults[assignment].push_back(elementIndex);
  }

  for (const auto& parsedResult : parsedResults) {
    if (std::find(std::begin(expectedResults), std::end(expectedResults),
                  parsedResult) == std::end(expectedResults))
      return -1;
  }
}

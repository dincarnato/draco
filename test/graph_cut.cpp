#include "graph_cut.hpp"

#include <armadillo>
#include <cassert>
#include <cmath>

void test_pairwise_distance() {
  arma::mat a{{1., 2.}, {3., 4.}, {5., 6.}};
  arma::mat b{{8., 12.}, {5., 21.}, {2., 3.}};

  auto distance =
      pairwise_distances(a, b.submat(arma::span::all, arma::span::all));
  arma::mat expected{{std::sqrt(149.), std::sqrt(377.), std::sqrt(2.)},
                     {std::sqrt(89.), std::sqrt(293.), std::sqrt(2.)},
                     {std::sqrt(45.), 15., std::sqrt(18.)}};
  assert(arma::approx_equal(distance, expected, "absdiff", 1e-12));
}

int main() { test_pairwise_distance(); }

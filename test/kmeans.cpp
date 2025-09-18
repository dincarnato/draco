#include "kmeans.hpp"

#include <armadillo>
#include <cassert>
#include <random>

void test_distance_squared_functions() {
  assert(kmeans::distance_squared(7.5f, 3.f) == 20.25f);
  assert(kmeans::distance_squared(7.5, 3.) == 20.25);

  arma::vec a{1., 2.5, 3.};
  arma::vec b{4.2, 0.7, -2.1};
  assert(kmeans::distance_squared(a, b) == 39.49);

  arma::mat m{{1., 2.5, 3.}, {4.2, 0.7, -2.1}};
  assert(kmeans::distance_squared(m.row(0), m.row(1)) == 39.49);
}

void test_kmeans_pp_centroids_initialization() {
  arma::mat points{
      {9., 9.5}, {0.7, 1.7}, {21., 2.2}, {1.5, 2.2},
      {20., 2.}, {1., 2.},   {8., 9.2},  {8.5, 10.},

  };

  std::mt19937 random_generator(42);
  auto centroids = kmeans::kmeans_pp_centroids_initialization(
      points.submat(arma::span::all, arma::span::all), 3, random_generator);
  assert(std::size(centroids) == 3);
  assert(centroids[0] == 1); // First point after the random shuffle
  assert(centroids[1] == 2);
  assert(centroids[2] == 0);
}

int main() {
  test_distance_squared_functions();
  test_kmeans_pp_centroids_initialization();
}

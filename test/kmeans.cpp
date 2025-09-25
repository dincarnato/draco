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

void test_kmeans() {
  arma::mat points{
      {9., 9.5}, {0.7, 1.7}, {21., 2.2}, {1.5, 2.2},
      {20., 2.}, {1., 2.},   {8., 9.2},  {8.5, 10.},

  };

  std::mt19937 random_generator(42);
  auto results = kmeans::run(points.submat(arma::span::all, arma::span::all), 3,
                             1, random_generator);
  assert(std::size(results.clusters_indices) == 3);
  {
    bool condition =
        results.clusters_indices[0] == std::vector{static_cast<arma::uword>(1),
                                                   static_cast<arma::uword>(3),
                                                   static_cast<arma::uword>(5)};
    assert(condition);
  }
  {
    bool condition =
        results.clusters_indices[1] ==
        std::vector{static_cast<arma::uword>(2), static_cast<arma::uword>(4)};
    assert(condition);
  }
  {
    bool condition =
        results.clusters_indices[2] == std::vector{static_cast<arma::uword>(0),
                                                   static_cast<arma::uword>(6),
                                                   static_cast<arma::uword>(7)};
    assert(condition);
  }

  {
    bool condition =
        arma::all(results.centroids.row(0) ==
                  arma::rowvec{(0.7 + 1.5 + 1.) / 3., (1.7 + 2.2 + 2.) / 3.});
    assert(condition);
  }
  {
    bool condition =
        arma::all(results.centroids.row(1) == arma::rowvec{20.5, 2.1});
    assert(condition);
  }
  {
    bool condition =
        arma::all(results.centroids.row(2) ==
                  arma::rowvec{(9. + 8. + 8.5) / 3., (9.5 + 9.2 + 10.) / 3.});
    assert(condition);
  }
}

void test_kmeans_trivial() {
  arma::mat points{
      {9., 9.5}, {0.7, 1.7}, {21., 2.2}, {1.5, 2.2},
      {20., 2.}, {1., 2.},   {8., 9.2},  {8.5, 10.},

  };

  arma::rowvec expected_centroid =
      arma::sum(points, 0) / static_cast<double>(points.n_rows);
  std::vector<arma::uword> expected_indices =
      std::views::iota(static_cast<arma::uword>(0), points.n_rows) |
      std::ranges::to<std::vector>();

  {
    std::mt19937 random_generator(42);
    auto results = kmeans::run(points.submat(arma::span::all, arma::span::all),
                               1, 1, random_generator);
    assert(results.centroids.n_rows == 1);
    assert(arma::all(results.centroids.row(0) == expected_centroid));
    assert(std::size(results.clusters_indices) == 1);
    assert(results.clusters_indices[0] == expected_indices);
  }

  {
    std::mt19937 random_generator(42);
    auto results = kmeans::run(points.submat(arma::span::all, arma::span::all),
                               0, 0, random_generator);
    assert(results.centroids.n_rows == 1);
    assert(arma::all(results.centroids.row(0) == expected_centroid));
    assert(std::size(results.clusters_indices) == 1);
    assert(results.clusters_indices[0] == expected_indices);
  }
}

void test_kmeans_multiple_iterations() {
  constexpr std::uint64_t random_seed = 1;
  arma::mat points{
      {9., 9.5}, {0.7, 1.7}, {21., 2.2}, {1.5, 2.2},
      {20., 2.}, {1., 2.},   {8., 9.2},  {8.5, 10.},

  };

  // Pre-check
  {
    std::mt19937 random_generator(random_seed);
    auto results = kmeans::run(points.submat(arma::span::all, arma::span::all),
                               3, 1, random_generator);
    assert(std::size(results.clusters_indices) == 3);
    {
      bool condition =
          results.clusters_indices[0] ==
          std::vector{static_cast<arma::uword>(0), static_cast<arma::uword>(1),
                      static_cast<arma::uword>(3), static_cast<arma::uword>(5),
                      static_cast<arma::uword>(6), static_cast<arma::uword>(7)};
      assert(condition);
    }
    {
      bool condition = results.clusters_indices[1] ==
                       std::vector{static_cast<arma::uword>(2)};
      assert(condition);
    }
    {
      bool condition =
          results.clusters_indices[2] == std::vector{
                                             static_cast<arma::uword>(4),
                                         };
      assert(condition);
    }
  }

  std::mt19937 random_generator(random_seed);
  auto results = kmeans::run(points.submat(arma::span::all, arma::span::all), 3,
                             10, random_generator);
  assert(std::size(results.clusters_indices) == 3);
  {
    bool condition =
        results.clusters_indices[0] ==
        std::vector{static_cast<arma::uword>(2), static_cast<arma::uword>(4)};
    assert(condition);
  }
  {
    bool condition =
        results.clusters_indices[1] == std::vector{static_cast<arma::uword>(1),
                                                   static_cast<arma::uword>(3),
                                                   static_cast<arma::uword>(5)};
    assert(condition);
  }
  {
    bool condition =
        results.clusters_indices[2] == std::vector{static_cast<arma::uword>(0),
                                                   static_cast<arma::uword>(6),
                                                   static_cast<arma::uword>(7)};
    assert(condition);
  }
}

int main() {
  test_distance_squared_functions();
  test_kmeans_pp_centroids_initialization();
  test_kmeans();
  test_kmeans_trivial();
  test_kmeans_multiple_iterations();
}

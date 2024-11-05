#include "ringmap_data.hpp"
#include "args.hpp"
#include "ringmap_matrix.hpp"
#include "ringmap_matrix_row.hpp"

#include <algorithm>
#include <array>
#include <initializer_list>
#include <iostream>
#include <ranges>

RingmapMatrixRow row(unsigned begin_index, unsigned end_index,
                     std::initializer_list<unsigned> modifiedIndices) {
  RingmapMatrixRow row{std::move(modifiedIndices)};
  row.unsafe_set_original_begin_index(begin_index);
  row.unsafe_set_original_end_index(end_index);
  return row;
}

void test_fractions_reads_by_reads_with_minimum_read_overlap() {
  RingmapMatrix matrix{11};
  matrix.addModifiedIndicesRow(row(2, 11, {2, 5, 7, 10}));
  matrix.addModifiedIndicesRow(row(3, 9, {4, 8}));
  matrix.addModifiedIndicesRow(row(1, 11, {3, 10}));
  matrix.addModifiedIndicesRow(row(1, 10, {1, 4, 9}));
  matrix.addModifiedIndicesRow(row(3, 9, {6, 7}));

  Args args;
  RingmapData ringmap_data{"AAAAAAAAAAA", std::move(matrix), 0, 11, args};
  auto window = ringmap_data.get_new_range(3, 9);
  assert(window.data().rows_size() == 5);
  WeightedClusters weights{6, 2};
  {
    auto &&cluster = weights.cluster(0);
    std::ranges::fill_n(std::ranges::begin(cluster), 3, 1.);
  }
  {
    auto &&cluster = weights.cluster(1);
    std::ranges::fill(cluster | std::views::drop(3), 1.);
  }

  {
    // Be sure about the behavior without minimum overlap required
    auto [clusters_fraction, clusters_pattern, clusters_assignment] =
        window.fractionReadsByWeights(weights, 0., 0);

    assert(clusters_assignment.size() == 5);
    std::ranges::all_of(
        std::views::iota(static_cast<unsigned int>(3)) |
            std::views::transform([&](auto row_index) {
              return clusters_assignment[window.data().getIndices(row_index)];
            }),
        [&](auto &&row) {
          return std::ranges::equal(row, std::array<unsigned short, 2>{1, 0});
        });
    std::ranges::all_of(
        std::views::iota(static_cast<unsigned int>(3),
                         static_cast<unsigned int>(6)) |
            std::views::transform([&](auto row_index) {
              return clusters_assignment[window.data().getIndices(row_index)];
            }),
        [&](auto &&row) {
          return std::ranges::equal(row, std::array<unsigned short, 2>{0, 1});
        });
  }

  auto [clusters_fraction, clusters_pattern, clusters_assignment] =
      window.fractionReadsByWeights(weights, 0.5, 0);

  assert(clusters_assignment.size() == 4);
  std::ranges::all_of(
      std::views::iota(static_cast<unsigned int>(2)) |
          std::views::transform([&](auto row_index) {
            return clusters_assignment[window.data().getIndices(row_index)];
          }),
      [&](auto &&row) {
        return std::ranges::equal(row, std::array<unsigned short, 2>{1, 0});
      });
  std::ranges::all_of(
      std::views::iota(static_cast<unsigned int>(2),
                       static_cast<unsigned int>(4)) |
          std::views::transform([&](auto row_index) {
            return clusters_assignment[window.data().getIndices(row_index)];
          }),
      [&](auto &&row) {
        return std::ranges::equal(row, std::array<unsigned short, 2>{0, 1});
      });
}

int main() {
  test_fractions_reads_by_reads_with_minimum_read_overlap();
  return 0;
}

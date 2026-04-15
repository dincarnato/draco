#include "compact_ringmap.hpp"
#include "mutation_map_transcript_read.hpp"
#include "ringmap_matrix.hpp"
#include "ringmap_matrix_row.hpp"
#include <algorithm>
#include <array>
#include <compare>
#include <iterator>
#include <ranges>
#include <string_view>

static std::array<MutationMapTranscriptRead, 12> reads{
    MutationMapTranscriptRead{.begin = 0u, .end = 18u, .indices = {1, 3, 16}},
    MutationMapTranscriptRead{.begin = 0u, .end = 18u, .indices = {}},
    MutationMapTranscriptRead{.begin = 1u, .end = 18u, .indices = {1, 3, 16}},
    MutationMapTranscriptRead{.begin = 5u, .end = 20u, .indices = {}},
    MutationMapTranscriptRead{.begin = 10u, .end = 28u, .indices = {20, 23}},
    MutationMapTranscriptRead{.begin = 0u, .end = 18u, .indices = {}},
    MutationMapTranscriptRead{.begin = 5u, .end = 18u, .indices = {6}},
    MutationMapTranscriptRead{.begin = 12u, .end = 28u, .indices = {20, 23}},
    MutationMapTranscriptRead{.begin = 14u, .end = 30u, .indices = {20, 23}},
    MutationMapTranscriptRead{.begin = 14u, .end = 30u, .indices = {20}},
    MutationMapTranscriptRead{
        .begin = 14u, .end = 30u, .indices = {20, 21, 23}},
    MutationMapTranscriptRead{.begin = 14u, .end = 30u, .indices = {19, 24}},
};

static void test_ringmap_matrix_row_helper_spaceship_operator_no_begin_end() {
  auto matrix_rows = reads | std::views::transform([](const auto &read) {
                       RingmapMatrixRow row;
                       std::ranges::copy(read.indices, std::back_inserter(row));

                       row.copy_begin_end_indices(MutationMapTranscriptRead{
                           .begin = 0u, .end = 30u, .indices = {}});
                       return row;
                     }) |
                     std::ranges::to<std::vector>();

  assert(RingmapMatrixRowHelper(matrix_rows[0]) <=>
             RingmapMatrixRowHelper(matrix_rows[2]) ==
         std::weak_ordering::equivalent);
  assert(RingmapMatrixRowHelper(matrix_rows[1]) <=>
             RingmapMatrixRowHelper(matrix_rows[3]) ==
         std::weak_ordering::equivalent);
  assert(RingmapMatrixRowHelper(matrix_rows[1]) <
         RingmapMatrixRowHelper(matrix_rows[0]));
  assert(RingmapMatrixRowHelper(matrix_rows[9]) <
         RingmapMatrixRowHelper(matrix_rows[10]));
  assert(RingmapMatrixRowHelper(matrix_rows[9]) <
         RingmapMatrixRowHelper(matrix_rows[8]));
  assert(RingmapMatrixRowHelper(matrix_rows[11]) <
         RingmapMatrixRowHelper(matrix_rows[8]));
}

static void test_ringmap_matrix_row_helper_spaceship_operator() {
  RingmapMatrixRow read_a;
  RingmapMatrixRow read_b;

  read_a.copy_begin_end_indices(
      MutationMapTranscriptRead{.begin = 5u, .end = 10u, .indices = {}});
  read_b.push_back(0u);
  read_b.copy_begin_end_indices(
      MutationMapTranscriptRead{.begin = 3u, .end = 10u, .indices = {}});
  assert(RingmapMatrixRowHelper(read_b) < RingmapMatrixRowHelper(read_a));

  read_b.clear();
  read_b.copy_begin_end_indices(
      MutationMapTranscriptRead{.begin = 5u, .end = 10u, .indices = {}});
  assert(RingmapMatrixRowHelper(read_b) <=> RingmapMatrixRowHelper(read_a) ==
         std::weak_ordering::equivalent);

  read_b.push_back(0u);
  read_b.copy_begin_end_indices(
      MutationMapTranscriptRead{.begin = 5u, .end = 7u, .indices = {}});
  assert(RingmapMatrixRowHelper(read_b) < RingmapMatrixRowHelper(read_a));

  read_b.copy_begin_end_indices(
      MutationMapTranscriptRead{.begin = 5u, .end = 10u, .indices = {}});
  assert(RingmapMatrixRowHelper(read_a) < RingmapMatrixRowHelper(read_b));
}

static void test_compact_ringmap_constructor() {
  constexpr std::string_view sequence = "AGCTAATTCCGCCGATTTATATGGACCATA";
  RingmapMatrix matrix(static_cast<std::uint32_t>(std::size(reads)),
                       std::size(sequence));
  for (auto const &read : reads) {
    matrix.addRead(read);
  }

  auto ringmap = CompactRingmap(matrix);
  assert(ringmap.max_modifications() == 3);
  assert(ringmap.n_rows() == 11);

  auto start_end_counts_sizes_and_modifications =
      ringmap.raw_start_end_count_sizes_and_modifications();
  std::array<MutationMapTranscriptRead const *, 11> expected_reads{
      &reads[1], &reads[0], &reads[2],  &reads[6], &reads[3],  &reads[4],
      &reads[7], &reads[9], &reads[11], &reads[8], &reads[10],
  };

  assert(std::ranges::is_sorted(expected_reads, [](auto a, auto b) {
    auto ordering = ([&] {
      auto ordering = a->begin <=> b->begin;
      if (ordering != std::weak_ordering::equivalent) {
        return ordering;
      }

      ordering = a->end <=> b->end;
      if (ordering != std::weak_ordering::equivalent) {
        return ordering;
      }

      ordering = std::size(a->indices) <=> std::size(b->indices);
      if (ordering != std::weak_ordering::equivalent) {
        return ordering;
      }

      return a->indices <=> b->indices;
    })();

    return ordering == std::weak_ordering::less;
  }));
  std::array<std::vector<unsigned> const *, 11> expected_modifications;
  std::ranges::copy(expected_reads |
                        std::views::transform(
                            [&](auto const &read) { return &read->indices; }),
                    std::ranges::begin(expected_modifications));

  for (auto [start_end_count_size_and_modifications, expected_indices,
             row_index] :
       std::views::zip(start_end_counts_sizes_and_modifications |
                           std::views::chunk(ringmap.row_size()),
                       expected_modifications,
                       std::views::iota(static_cast<std::uint32_t>(0)))) {
    auto row = ringmap.row(row_index);

    auto expected_count = std::ranges::count_if(reads, [&](auto &read) {
      return read.begin == start_end_count_size_and_modifications[0] and
             read.end == start_end_count_size_and_modifications[1] and
             std::ranges::equal(read.indices, *expected_indices);
    });
    assert(start_end_count_size_and_modifications[0] == row.begin_index());
    assert(start_end_count_size_and_modifications[1] == row.end_index());
    assert(start_end_count_size_and_modifications[2] == expected_count);
    assert(row.count() == expected_count);
    assert(start_end_count_size_and_modifications[3] ==
           std::size(*expected_indices));
    assert(row.size() == std::size(*expected_indices));
    assert(std::ranges::equal(
        start_end_count_size_and_modifications | std::views::drop(4) |
            std::views::take(start_end_count_size_and_modifications[3]),
        *expected_indices));
    assert(std::ranges::equal(row.indices(), *expected_indices));
  }

  std::array<std::vector<std::uint32_t>, 11> expected_mapping{
      std::vector{1u, 5u}, std::vector{0u}, std::vector{2u},  std::vector{6u},
      std::vector{3u},     std::vector{4u}, std::vector{7u},  std::vector{9u},
      std::vector{11u},    std::vector{8u}, std::vector{10u},
  };
  assert(std::ranges::equal(
      ringmap.mapping(), expected_mapping,
      [](auto const &a, auto const &b) { return std::ranges::equal(a, b); }));
}

static void test_compact_ringmap_iterator() {
  constexpr std::string_view sequence = "AGCTAATTCCGCCGATTTATATGGACCATA";
  RingmapMatrix matrix(static_cast<std::uint32_t>(std::size(reads)),
                       std::size(sequence));
  for (auto const &read : reads) {
    matrix.addRead(read);
  }

  auto ringmap = CompactRingmap(matrix);
  for (auto [row_index, row] : std::views::zip(
           std::views::iota(static_cast<std::uint32_t>(0)), ringmap)) {
    assert(row == ringmap.row(row_index));
  }

  assert(std::ranges::begin(ringmap)[4] == ringmap.row(4));
  assert(*(std::ranges::end(ringmap) - 2) == ringmap.row(ringmap.n_rows() - 2));
}

int main() {
  test_ringmap_matrix_row_helper_spaceship_operator_no_begin_end();
  test_ringmap_matrix_row_helper_spaceship_operator();
  test_compact_ringmap_constructor();
  test_compact_ringmap_iterator();
}

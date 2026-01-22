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

static void test_ringmap_matrix_row_helper_spaceship_operator() {
  auto matrix_rows = reads | std::views::transform([](const auto &read) {
                       RingmapMatrixRow row;
                       std::ranges::copy(read.indices, std::back_inserter(row));
                       row.copy_begin_end_indices(read);
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

static void test_compact_ringmap_constructor() {
  constexpr std::string_view sequence = "AGCTAATTCCGCCGATTTATATGGACCATA";
  RingmapMatrix matrix(static_cast<std::uint32_t>(std::size(reads)),
                       std::size(sequence));
  for (auto const &read : reads) {
    matrix.addRead(read);
  }

  auto ringmap = CompactRingmap(matrix);
  assert(ringmap.max_modifications() == 3);
  assert(ringmap.n_rows() == 7);

  auto counts_sizes_and_modifications =
      ringmap.raw_count_sizes_and_modifications();
  std::array<std::vector<unsigned> const *, 7> expected_modifications{
      &reads[1].indices,  &reads[6].indices, &reads[9].indices,
      &reads[11].indices, &reads[4].indices, &reads[0].indices,
      &reads[10].indices,
  };
  assert(std::ranges::is_sorted(expected_modifications, [](auto a, auto b) {
    auto size_a = std::size(*a);
    auto size_b = std::size(*b);
    if (size_a == size_b) {
      return *a < *b;
    } else {
      return size_a < size_b;
    }
  }));
  for (auto [count_size_and_modifications, expected_indices, row_index] :
       std::views::zip(counts_sizes_and_modifications |
                           std::views::chunk(ringmap.row_size()),
                       expected_modifications,
                       std::views::iota(static_cast<std::uint32_t>(0)))) {

    auto expected_count = std::ranges::count_if(reads, [&](auto &read) {
      return std::ranges::equal(read.indices, *expected_indices);
    });
    assert(count_size_and_modifications[0] == expected_count);
    assert(count_size_and_modifications[1] == std::size(*expected_indices));
    assert(std::ranges::equal(
        count_size_and_modifications | std::views::drop(2) |
            std::views::take(count_size_and_modifications[1]),
        *expected_indices));
  }
}

int main() {
  test_ringmap_matrix_row_helper_spaceship_operator();
  test_compact_ringmap_constructor();
}

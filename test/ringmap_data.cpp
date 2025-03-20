#include "ringmap_data.hpp"
#include "ringmap_matrix.hpp"
#include <algorithm>
#include <iterator>
#include <ranges>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

RingmapMatrix generate_ringmap_matrix(std::string const &sequence,
                                      std::span<const unsigned> low_freq_bases,
                                      unsigned n_reads, unsigned n_bases,
                                      unsigned single_mutation_reads,
                                      unsigned double_mutation_reads) {
  auto modificable_bases =
      std::views::zip(std::views::iota(0u), sequence) |
      std::views::filter([](auto pair) {
        auto base = std::get<1>(pair);
        return base == 'A' or base == 'C';
      }) |
      std::views::transform([](auto pair) { return std::get<0>(pair); }) |
      std::ranges::to<std::vector>();
  auto n_modificable_bases =
      static_cast<unsigned>(std::ranges::distance(modificable_bases));

  auto filter_base = [&](unsigned iter_index, unsigned base_index) {
    return not std::ranges::binary_search(low_freq_bases, base_index) or
           iter_index % 6 == 0;
  };

  RingmapMatrix matrix(n_reads, n_bases);
  for (auto iter_index : std::views::iota(0u, single_mutation_reads)) {
    std::vector<unsigned> indices;
    if (filter_base(iter_index, iter_index % n_bases)) {
      indices.emplace_back(modificable_bases[iter_index % n_modificable_bases]);
    }
    matrix.addRead(MutationMapTranscriptRead{
        .begin = 0u,
        .end = n_bases,
        .indices = std::move(indices),
    });
  }
  for (auto iter_index : std::views::iota(0u, double_mutation_reads)) {
    unsigned base_a = modificable_bases[(iter_index / n_modificable_bases) %
                                        n_modificable_bases];
    unsigned base_b = modificable_bases[iter_index % n_modificable_bases];

    if (base_b < base_a) {
      std::swap(base_a, base_b);
    }

    std::vector<unsigned> indices;
    if (filter_base(iter_index, base_a)) {
      indices.emplace_back(base_a);
    }
    if (filter_base(iter_index, base_b)) {
      indices.emplace_back(base_b);
    }

    matrix.addRead(MutationMapTranscriptRead{
        .begin = 0u,
        .end = n_bases,
        .indices = std::move(indices),
    });
  }

  for ([[maybe_unused]] auto iter_index : std::views::iota(
           single_mutation_reads + double_mutation_reads, n_reads)) {
    matrix.addRead(MutationMapTranscriptRead{
        .begin = 0u,
        .end = n_bases,
        .indices = std::vector<unsigned>(),
    });
  }

  return matrix;
}

void test_filter_bases_on_replicates_same_data() {
  constexpr std::string_view sequence_view = "ACATATTACGGAGGGTACTT";
  constexpr auto n_reads = 10000u;
  constexpr auto n_bases = static_cast<unsigned>(std::size(sequence_view));
  constexpr auto low_freq_bases = std::array{2u, 8u, 16u};

  std::string sequence(sequence_view);
  auto matrix = generate_ringmap_matrix(sequence, low_freq_bases, n_reads,
                                        n_bases, 100u, 900u);

  Args args;
  auto ringmaps_data =
      std::views::iota(0, 3) | std::views::transform([&](auto) {
        return RingmapData(sequence, RingmapMatrix(matrix), 0, n_bases, args);
      }) |
      std::ranges::to<std::vector<RingmapData>>();

  auto single_filtered_ringmap = RingmapData(ringmaps_data[0]);
  single_filtered_ringmap.filterBases();
  assert(single_filtered_ringmap != ringmaps_data[0]);

  RingmapData::filter_bases_on_replicates(ringmaps_data);

  assert(std::ranges::all_of(ringmaps_data, [&](const auto &ringmap_data) {
    return single_filtered_ringmap == ringmap_data;
  }));
}

void test_filter_bases_on_replicates_same_data_refilter() {
  constexpr std::string_view sequence_view = "ACATATTACGGAGGGTACTT";
  constexpr auto n_reads = 10000u;
  constexpr auto n_bases = static_cast<unsigned>(std::size(sequence_view));
  constexpr auto low_freq_bases = std::array{2u, 8u, 16u};

  std::string sequence(sequence_view);

  auto remove_base_5_mutations = [](RingmapData &ringmap_data) {
    auto &data = test::RingmapData::m_data(ringmap_data);
    auto rows = data.rows();
    for (auto &&data_row : rows) {
      data_row[2] = false;
    }
  };

  auto matrix = generate_ringmap_matrix(sequence, low_freq_bases, n_reads,
                                        n_bases, 100u, 900u);
  Args args;

  auto ringmaps_data =
      std::views::iota(0, 3) | std::views::transform([&](auto) {
        return RingmapData(sequence, RingmapMatrix(matrix), 0, n_bases, args);
      }) |
      std::ranges::to<std::vector<RingmapData>>();

  auto single_filtered_ringmap = RingmapData(ringmaps_data[0]);
  single_filtered_ringmap.filterBases();
  RingmapData::filter_bases_on_replicates(ringmaps_data);

  remove_base_5_mutations(single_filtered_ringmap);
  for (auto &ringmap_data : ringmaps_data) {
    remove_base_5_mutations(ringmap_data);
  }

  single_filtered_ringmap.filterBases();
  RingmapData::filter_bases_on_replicates(ringmaps_data);

  assert(std::ranges::all_of(ringmaps_data, [&](const auto &ringmap_data) {
    return single_filtered_ringmap == ringmap_data;
  }));
}

void test_filter_bases_on_replicates_different_data() {
  constexpr std::string_view sequence_view = "ACATATTACGGAGGGTACTT";
  constexpr auto n_reads = 10000u;
  constexpr auto n_bases = static_cast<unsigned>(std::size(sequence_view));

  std::string sequence(sequence_view);
  std::vector<RingmapData> ringmaps_data;

  Args args;
  {
    constexpr auto low_freq_bases = std::array{2u, 8u, 16u};
    auto matrix = generate_ringmap_matrix(sequence, low_freq_bases, n_reads,
                                          n_bases, 100u, 900u);

    ringmaps_data.emplace_back(sequence, RingmapMatrix(matrix), 0, n_bases,
                               args);
  }
  {
    constexpr auto low_freq_bases = std::array{1u, 2u, 9u};
    auto matrix = generate_ringmap_matrix(sequence, low_freq_bases, n_reads,
                                          n_bases, 100u, 900u);

    ringmaps_data.emplace_back(sequence, RingmapMatrix(matrix), 0, n_bases,
                               args);
  }

  auto singly_filtered_ringmaps =
      ringmaps_data |
      std::views::transform([&](auto const &original_ringmap_data) {
        auto ringmap_data = RingmapData(original_ringmap_data);
        ringmap_data.filterBases();
        return ringmap_data;
      }) |
      std::ranges::to<std::vector>();
  assert(singly_filtered_ringmaps[0] != singly_filtered_ringmaps[1]);

  RingmapData::filter_bases_on_replicates(ringmaps_data);

  assert(singly_filtered_ringmaps[0] != ringmaps_data[0]);
  assert(singly_filtered_ringmaps[0] != ringmaps_data[1]);
  assert(singly_filtered_ringmaps[1] != ringmaps_data[0]);
  assert(singly_filtered_ringmaps[1] != ringmaps_data[1]);

  assert(
      singly_filtered_ringmaps[0].getNonFilteredToFilteredMap().contains(1u));
  assert(not singly_filtered_ringmaps[0].getNonFilteredToFilteredMap().contains(
      2u));

  constexpr auto all_low_freq_bases = std::array{1u, 2u, 8u, 9u, 16u};
  assert(std::ranges::all_of(ringmaps_data, [&](const auto &ringmap_data) {
    auto const &non_filtered_to_filtered =
        ringmap_data.getNonFilteredToFilteredMap();
    return std::ranges::all_of(
        std::views::zip(std::views::iota(0u), sequence_view), [&](auto &&pair) {
          auto &&[base_index, base] = pair;
          auto should_be_there =
              (base == 'A' or base == 'C') and
              not std::ranges::contains(all_low_freq_bases, base_index);
          return non_filtered_to_filtered.contains(base_index) ==
                 should_be_there;
        });
  }));
}

int main() {
  test_filter_bases_on_replicates_same_data();
  test_filter_bases_on_replicates_same_data_refilter();
  test_filter_bases_on_replicates_different_data();
}

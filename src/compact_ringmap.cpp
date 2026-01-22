#include "compact_ringmap.hpp"
#include "ringmap_matrix.hpp"
#include <algorithm>
#include <compare>
#include <iterator>
#include <map>
#include <memory>
#include <ranges>
#include <span>
#include <vector>

CompactRingmap::CompactRingmap(RingmapMatrix const &ringmap_matrix) {
  std::map<RingmapMatrixRowHelper, std::uint32_t> unique_modifications_sets;
  std::uint32_t max_modifications = 0;

  for (auto rows = ringmap_matrix.rows(); auto &&row : rows) {
    auto modified_indices = RingmapMatrixRowHelper(row.modifiedIndices());
    max_modifications = std::max(
        max_modifications,
        static_cast<std::uint32_t>(std::size(*modified_indices.inner)));

    if (auto iter = unique_modifications_sets.find(modified_indices);
        iter == std::ranges::end(unique_modifications_sets)) {
      unique_modifications_sets.insert(std::pair{modified_indices, 1});
    } else {
      iter->second += 1;
    }
  }

  auto n_rows =
      static_cast<std::uint32_t>(std::size(unique_modifications_sets));
  auto row_size = max_modifications + 2;
  auto storage_size = n_rows * row_size;
  std::unique_ptr<std::uint32_t[]> sizes_and_modifications(
      new std::uint32_t[storage_size]);
  std::ranges::for_each(
      std::views::zip(unique_modifications_sets,
                      std::span(sizes_and_modifications.get(), storage_size) |
                          std::views::chunk(row_size)),
      [&](auto &&tuple) {
        auto [row_with_count, size_and_modifications] = tuple;
        auto [row, count] = row_with_count;

        size_and_modifications[0] = count;
        size_and_modifications[1] =
            static_cast<std::uint32_t>(std::size(*row.inner));
        std::ranges::copy(
            *row.inner | std::views::transform(
                             [&](auto base_index) { return base_index; }),
            std::next(std::ranges::begin(size_and_modifications), 2));
      });

  count_sizes_and_modifications_ = std::move(sizes_and_modifications);
  n_rows_ = n_rows;
  max_modifications_ = max_modifications;
}

RingmapMatrixRowHelper::RingmapMatrixRowHelper(
    RingmapMatrixRow const &row) noexcept
    : inner(&row) {}

std::weak_ordering RingmapMatrixRowHelper::operator<=>(
    RingmapMatrixRowHelper const &other) const noexcept {
  auto size_ordering = std::size(*inner) <=> std::size(*other.inner);
  if (size_ordering != std::weak_ordering::equivalent) {
    return size_ordering;
  }

  for (auto [a, b] : std::views::zip(*inner, *other.inner)) {
    auto cmp = a <=> b;
    if (cmp != std::weak_ordering::equivalent) {
      return cmp;
    }
  }
  return std::weak_ordering::equivalent;
}

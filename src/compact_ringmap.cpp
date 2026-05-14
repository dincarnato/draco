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

CompactRingmap::CompactRingmap(RingmapMatrix const &ringmap_matrix,
                               std::uint32_t start_index,
                               std::uint32_t end_index)
    : original_n_rows_(ringmap_matrix.rows_size()), start_index_(start_index),
      end_index_(end_index) {
  std::map<RingmapMatrixRowHelper,
           std::tuple<std::uint32_t, std::vector<std::uint32_t>>>
      unique_modifications_sets;
  std::uint32_t max_modifications = 0;

  for (auto rows =
           std::views::zip(std::views::iota(static_cast<std::uint32_t>(0)),
                           ringmap_matrix.rows());
       auto &&[row_index, row] : rows) {
    auto modified_indices = RingmapMatrixRowHelper(row.modifiedIndices());
    max_modifications = std::max(
        max_modifications,
        static_cast<std::uint32_t>(std::size(*modified_indices.inner)));

    if (auto iter = unique_modifications_sets.find(modified_indices);
        iter == std::ranges::end(unique_modifications_sets)) {
      unique_modifications_sets.insert(
          std::pair{modified_indices, std::tuple{1, std::vector{row_index}}});
    } else {
      auto &&stored_data = iter->second;
      std::get<0>(stored_data) += 1;
      std::get<1>(stored_data).push_back(row_index);
    }
  }

  auto n_rows =
      static_cast<std::uint32_t>(std::size(unique_modifications_sets));
  auto row_size = max_modifications + 4;
  auto storage_size = n_rows * row_size;
  std::unique_ptr<std::uint32_t[]> start_end_count_sizes_and_modifications(
      new std::uint32_t[storage_size]);
  std::vector<std::vector<std::uint32_t>> mapping;
  mapping.reserve(n_rows);
  std::ranges::for_each(
      std::views::zip(unique_modifications_sets,
                      std::span(start_end_count_sizes_and_modifications.get(),
                                storage_size) |
                          std::views::chunk(row_size)),
      [&](auto &&tuple) {
        auto [row_with_count_and_indices,
              start_end_count_size_and_modifications] = tuple;
        auto [row, count_and_indices] = row_with_count_and_indices;
        auto [count, indices] = count_and_indices;

        start_end_count_size_and_modifications[0] = row.inner->begin_index();
        start_end_count_size_and_modifications[1] = row.inner->end_index();
        start_end_count_size_and_modifications[2] = count;
        start_end_count_size_and_modifications[3] =
            static_cast<std::uint32_t>(std::size(*row.inner));
        std::ranges::copy(
            *row.inner | std::views::transform(
                             [&](auto base_index) { return base_index; }),
            std::next(
                std::ranges::begin(start_end_count_size_and_modifications), 4));
        mapping.push_back(std::move(indices));
      });

  start_end_count_sizes_and_modifications_ =
      std::move(start_end_count_sizes_and_modifications);
  n_rows_ = n_rows;
  max_modifications_ = max_modifications;
  mapping_ = std::move(mapping);
}

RingmapMatrixRowHelper::RingmapMatrixRowHelper(
    RingmapMatrixRow const &row) noexcept
    : inner(&row) {}

std::weak_ordering RingmapMatrixRowHelper::operator<=>(
    RingmapMatrixRowHelper const &other) const noexcept {
  auto ordering = inner->begin_index() <=> other.inner->begin_index();
  if (ordering != std::weak_ordering::equivalent) {
    return ordering;
  }

  ordering = inner->end_index() <=> other.inner->end_index();
  if (ordering != std::weak_ordering::equivalent) {
    return ordering;
  }

  ordering = std::size(*inner) <=> std::size(*other.inner);
  if (ordering != std::weak_ordering::equivalent) {
    return ordering;
  }

  for (auto [a, b] : std::views::zip(*inner, *other.inner)) {
    auto cmp = a <=> b;
    if (cmp != std::weak_ordering::equivalent) {
      return cmp;
    }
  }
  return std::weak_ordering::equivalent;
}

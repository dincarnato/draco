#pragma once

#include <cassert>
#include <compare>
#include <cstdint>
#include <memory>
#include <span>

class RingmapMatrix;
struct RingmapMatrixRow;

/**
 * A compact version of the modifications matrix
 *
 * This abstraction does **not** track the begin and end of the single reads,
 * it only stores the modified indices. In this way it is possible to keep a
 * count for each set of modifications and improve the memory footprint.
 */
struct CompactRingmap {
  CompactRingmap() = default;
  explicit CompactRingmap(RingmapMatrix const &ringmap_matrix);

  constexpr std::uint32_t row_size() const noexcept {
    return max_modifications_ + 2;
  }

  constexpr std::uint32_t n_rows() const noexcept { return n_rows_; }

  constexpr std::uint32_t max_modifications() const noexcept {
    return max_modifications_;
  }

  constexpr std::span<std::uint32_t const>
  raw_count_sizes_and_modifications() const noexcept {
    return std::span(count_sizes_and_modifications_.get(),
                     row_size() * n_rows_);
  }

  constexpr std::span<std::uint32_t const>
  raw_row(std::uint32_t row) const noexcept {
    assert(row < n_rows_);
    auto begin = row * row_size();
    return std::span(count_sizes_and_modifications_.get() + begin, row_size());
  }

protected:
  std::unique_ptr<std::uint32_t[]> count_sizes_and_modifications_{};
  std::uint32_t n_rows_{};
  std::uint32_t max_modifications_{};
};

struct RingmapMatrixRowHelper {
  RingmapMatrixRowHelper(RingmapMatrixRow const &row) noexcept;

  std::weak_ordering
  operator<=>(RingmapMatrixRowHelper const &other) const noexcept;

  RingmapMatrixRow const *inner;
};

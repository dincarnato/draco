#pragma once

#include <cassert>
#include <compare>
#include <cstdint>
#include <iterator>
#include <memory>
#include <span>
#include <vector>

class RingmapMatrix;
struct RingmapMatrixRow;

struct CompactRingmap;

struct CompactRingmapRow {
  constexpr CompactRingmapRow(CompactRingmap const &compact_ringmap,
                              std::uint32_t row) noexcept
      : compact_ringmap_(&compact_ringmap), row_(row) {}

  constexpr std::span<std::uint32_t const> raw() const noexcept;
  constexpr std::span<std::uint32_t const> mapped_rows() const noexcept;

  constexpr std::uint32_t begin_index() const noexcept { return raw()[0]; }
  constexpr std::uint32_t end_index() const noexcept { return raw()[1]; }
  constexpr std::uint32_t read_size() const noexcept {
    return end_index() - begin_index();
  }
  constexpr std::uint32_t count() const noexcept { return raw()[2]; }
  constexpr std::uint32_t size() const noexcept { return raw()[3]; }
  constexpr std::span<std::uint32_t const> indices() const noexcept {
    return raw().subspan(4, size());
  }

  constexpr bool operator==(CompactRingmapRow const &other) const noexcept {
    return compact_ringmap_ == other.compact_ringmap_ && row_ == other.row_;
  }

protected:
  CompactRingmap const *compact_ringmap_;
  std::uint32_t row_;
};

struct CompactRingmapIterator {
  friend struct CompactRingmapRange;

  using value_type = CompactRingmapRow;
  using self = CompactRingmapIterator;

  constexpr CompactRingmapIterator() noexcept = default;
  constexpr CompactRingmapIterator(CompactRingmap const &compact_ringmap,
                                   std::uint32_t row) noexcept
      : compact_ringmap_(&compact_ringmap), row_(row) {}

  constexpr std::partial_ordering
  operator<=>(self const &other) const noexcept {
    if (compact_ringmap_ != other.compact_ringmap_) {
      return std::partial_ordering::unordered;
    }

    return row_ <=> other.row_;
  }
  constexpr bool operator==(self const &other) const noexcept = default;

  constexpr value_type operator*() const noexcept {
    return CompactRingmapRow(*compact_ringmap_, row_);
  }

  constexpr value_type operator[](std::int64_t offset) const noexcept {
    return CompactRingmapRow(
        *compact_ringmap_,
        static_cast<std::uint32_t>(static_cast<std::int64_t>(row_) + offset));
  }

  constexpr self &operator++() noexcept {
    ++row_;
    return *this;
  }

  constexpr self operator++(int) noexcept {
    auto other = *this;
    ++row_;
    return other;
  }

  constexpr self &operator--() noexcept {
    --row_;
    return *this;
  }

  constexpr self operator--(int) noexcept {
    auto other = *this;
    --row_;
    return other;
  }

  constexpr self &operator+=(std::int64_t offset) noexcept {
    row_ = static_cast<std::uint32_t>(static_cast<std::int64_t>(row_) + offset);
    return *this;
  }

  constexpr self &operator-=(std::int64_t offset) noexcept {
    row_ = static_cast<std::uint32_t>(static_cast<std::int64_t>(row_) - offset);
    return *this;
  }

  constexpr self operator+(std::int64_t offset) const noexcept {
    auto other = *this;
    other += offset;
    return other;
  }

  friend constexpr self operator+(std::int64_t offset,
                                  self const &it) noexcept {
    return it + offset;
  }

  constexpr self operator-(std::int64_t offset) const noexcept {
    auto other = *this;
    other -= offset;
    return other;
  }

  constexpr std::int64_t operator-(self const &other) const noexcept {
    return static_cast<std::int64_t>(row_) -
           static_cast<std::int64_t>(other.row_);
  }

protected:
  CompactRingmap const *compact_ringmap_{};
  std::uint32_t row_{};
};

static_assert(std::random_access_iterator<CompactRingmapIterator>);

struct CompactRingmapRange {
  using iterator = CompactRingmapIterator;

  constexpr CompactRingmapRange() noexcept = default;
  constexpr CompactRingmapRange(iterator begin, iterator end) noexcept
      : compact_ringmap_(begin.compact_ringmap_), begin_row_(begin.row_),
        end_row_(end.row_) {
    assert(begin.compact_ringmap_ == end.compact_ringmap_);
  }

  constexpr iterator begin() const noexcept {
    return CompactRingmapIterator(*compact_ringmap_, begin_row_);
  }

  constexpr iterator end() const noexcept {
    return CompactRingmapIterator(*compact_ringmap_, end_row_);
  }

protected:
  CompactRingmap const *compact_ringmap_{};
  std::uint32_t begin_row_{};
  std::uint32_t end_row_{};
};

/**
 * A compact version of the modifications matrix
 *
 * This abstraction does **not** track the begin and end of the single reads,
 * it only stores the modified indices. In this way it is possible to keep a
 * count for each set of modifications and improve the memory footprint.
 */
struct CompactRingmap {
  CompactRingmap() = default;
  explicit CompactRingmap(RingmapMatrix const &ringmap_matrix,
                          std::uint32_t start_index, std::uint32_t end_index);

  constexpr std::uint32_t row_size() const noexcept {
    return max_modifications_ + 4;
  }

  constexpr std::uint32_t n_rows() const noexcept { return n_rows_; }

  constexpr std::uint32_t max_modifications() const noexcept {
    return max_modifications_;
  }

  constexpr std::uint32_t original_n_rows() const noexcept {
    return original_n_rows_;
  }

  constexpr std::span<const std::vector<std::uint32_t>>
  mapping() const noexcept {
    return mapping_;
  }

  constexpr std::span<std::uint32_t const>
  raw_start_end_count_sizes_and_modifications() const noexcept {
    return std::span(start_end_count_sizes_and_modifications_.get(),
                     row_size() * n_rows_);
  }

  constexpr std::span<std::uint32_t const>
  raw_row(std::uint32_t row) const noexcept {
    assert(row < n_rows_);
    auto begin = row * row_size();
    return std::span(start_end_count_sizes_and_modifications_.get() + begin,
                     row_size());
  }

  constexpr CompactRingmapRow row(std::uint32_t row) const noexcept {
    return CompactRingmapRow(*this, row);
  }

  constexpr CompactRingmapIterator begin() const noexcept {
    return CompactRingmapIterator(*this, 0);
  }

  constexpr CompactRingmapIterator end() const noexcept {
    return CompactRingmapIterator(*this, n_rows_);
  }

  constexpr std::uint32_t start_index() const noexcept { return start_index_; }
  constexpr std::uint32_t end_index() const noexcept { return end_index_; }

protected:
  std::unique_ptr<std::uint32_t[]> start_end_count_sizes_and_modifications_{};
  std::vector<std::vector<std::uint32_t>> mapping_;
  std::uint32_t n_rows_{};
  std::uint32_t max_modifications_{};
  std::uint32_t original_n_rows_{};
  std::uint32_t start_index_{};
  std::uint32_t end_index_{};
};

struct RingmapMatrixRowHelper {
  RingmapMatrixRowHelper(RingmapMatrixRow const &row) noexcept;

  std::weak_ordering
  operator<=>(RingmapMatrixRowHelper const &other) const noexcept;

  RingmapMatrixRow const *inner;
};

constexpr std::span<std::uint32_t const>
CompactRingmapRow::raw() const noexcept {
  return compact_ringmap_->raw_row(row_);
}

constexpr std::span<std::uint32_t const>
CompactRingmapRow::mapped_rows() const noexcept {
  return compact_ringmap_->mapping()[row_];
}

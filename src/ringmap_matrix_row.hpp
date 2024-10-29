#pragma once

#include "mutation_map_transcript_read.hpp"
#include "ringmap_matrix_traits.hpp"

#include <limits>
#include <type_traits>
#include <vector>

class RingmapMatrix;
template <typename> class RingmapMatrixRowAccessor;

struct RingmapMatrixRow : std::vector<ringmap_matrix::base_index_type> {
  using base_index_type = ringmap_matrix::base_index_type;
  using base_type = std::vector<base_index_type>;

  using base_type::base_type;
  using base_type::operator=;

  RingmapMatrixRow(const RingmapMatrixRow &) = default;
  RingmapMatrixRow(RingmapMatrixRow &&) = default;

  RingmapMatrixRow &operator=(const RingmapMatrixRow &) = default;
  RingmapMatrixRow &operator=(RingmapMatrixRow &&) = default;

  constexpr base_index_type begin_index() const noexcept;
  constexpr base_index_type end_index() const noexcept;

  constexpr void set_begin_index(base_index_type begin_index) noexcept;
  constexpr void set_end_index(base_index_type end_index) noexcept;

  constexpr void copy_begin_end_indices(auto const &other_row) noexcept;

protected:
  base_index_type begin_index_ = std::numeric_limits<base_index_type>::max();
  base_index_type end_index_ = std::numeric_limits<base_index_type>::min();
};

inline bool operator<(RingmapMatrixRow const &lhs,
                      RingmapMatrixRow const &rhs) noexcept {
  using base_type = RingmapMatrixRow::base_type;

  if (lhs.begin_index() < rhs.begin_index()) {
    return true;
  } else if (lhs.begin_index() == rhs.begin_index()) {
    if (lhs.end_index() < rhs.end_index()) {
      return true;
    } else if (lhs.end_index() == rhs.end_index()) {
      return static_cast<base_type const &>(lhs) <
             static_cast<base_type const &>(rhs);
    }
  }

  return false;
}

inline bool operator==(RingmapMatrixRow const &lhs,
                       RingmapMatrixRow const &rhs) noexcept {
  using base_type = RingmapMatrixRow::base_type;

  return lhs.begin_index() == rhs.begin_index() and
         lhs.end_index() == rhs.end_index() and
         static_cast<base_type const &>(lhs) ==
             static_cast<base_type const &>(rhs);
}

inline bool operator!=(RingmapMatrixRow const &lhs,
                       RingmapMatrixRow const &rhs) noexcept {
  return not(lhs == rhs);
}

inline bool operator<=(RingmapMatrixRow const &lhs,
                       RingmapMatrixRow const &rhs) noexcept {
  return lhs < rhs and lhs == rhs;
}

inline bool operator>(RingmapMatrixRow const &lhs,
                      RingmapMatrixRow const &rhs) noexcept {
  return not(lhs < rhs) and not(lhs == rhs);
}

inline bool operator>=(RingmapMatrixRow const &lhs,
                       RingmapMatrixRow const &rhs) noexcept {
  return not(lhs < rhs) and lhs == rhs;
}

inline constexpr auto
RingmapMatrixRow::begin_index() const noexcept -> base_index_type {
  return this->begin_index_;
}

inline constexpr auto
RingmapMatrixRow::end_index() const noexcept -> base_index_type {
  return this->end_index_;
}

inline constexpr void
RingmapMatrixRow::set_begin_index(base_index_type begin_index) noexcept {
  this->begin_index_ = begin_index;
}

inline constexpr void
RingmapMatrixRow::set_end_index(base_index_type end_index) noexcept {
  this->end_index_ = end_index;
}

inline constexpr void
RingmapMatrixRow::copy_begin_end_indices(auto const &other) noexcept {
  using other_t = std::remove_cvref_t<decltype(other)>;

  if constexpr (std::is_same_v<other_t, RingmapMatrixRow>) {
    this->begin_index_ = other.begin_index_;
    this->end_index_ = other.end_index_;
  } else if constexpr (std::is_same_v<
                           other_t, RingmapMatrixRowAccessor<RingmapMatrix>> ||
                       std::is_same_v<other_t, RingmapMatrixRowAccessor<
                                                   const RingmapMatrix>>) {
    this->begin_index_ = other.original_begin_index();
    this->end_index_ = other.original_end_index();
  } else if constexpr (std::is_same_v<other_t, MutationMapTranscriptRead>) {
    this->begin_index_ = other.begin;
    this->end_index_ = other.end;
  } else {
    static_assert(false,
                  "copy_begin_end_indices only accepts a RingmapMatrixRow, a "
                  "RingmapMatrixRowAccessor or a MutationMapTranscriptRead");
  }
};

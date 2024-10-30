#pragma once

#include "mutation_map_transcript_read.hpp"
#include "ringmap_matrix_traits.hpp"

#include <concepts>
#include <limits>
#include <type_traits>
#include <vector>

class RingmapMatrix;
template <typename> class RingmapMatrixRowAccessor;

namespace results {
struct Window;
} /* namespace results */

struct RingmapMatrixWindowIndices {
  ringmap_matrix::base_index_type begin_index;
  ringmap_matrix::base_index_type end_index;
};

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

  constexpr base_index_type original_begin_index() const noexcept;
  constexpr base_index_type original_end_index() const noexcept;

  constexpr base_index_type window_begin_index() const noexcept;
  constexpr base_index_type window_end_index() const noexcept;

  constexpr void
  unsafe_set_original_begin_index(base_index_type begin_index) noexcept;
  constexpr void
  unsafe_set_original_end_index(base_index_type end_index) noexcept;

  constexpr void copy_begin_end_indices(auto const &other_row) noexcept;
  constexpr void copy_window_begin_end_indices(auto const &window) noexcept;

protected:
  base_index_type begin_index_ = std::numeric_limits<base_index_type>::max();
  base_index_type end_index_ = std::numeric_limits<base_index_type>::min();
  base_index_type window_begin_index_ = 0;
  base_index_type window_end_index_ =
      std::numeric_limits<base_index_type>::max();
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
  if (this->begin_index_ == std::numeric_limits<base_index_type>::max()) {
    return std::numeric_limits<base_index_type>::max();
  } else {
    return std::max(this->begin_index_, this->window_begin_index_);
  }
}

inline constexpr auto
RingmapMatrixRow::end_index() const noexcept -> base_index_type {
  if (this->end_index_ == std::numeric_limits<base_index_type>::min()) {
    return std::numeric_limits<base_index_type>::min();
  } else {
    return std::min(this->end_index_, this->window_end_index_);
  }
}

inline constexpr auto
RingmapMatrixRow::original_begin_index() const noexcept -> base_index_type {
  return this->begin_index_;
}

inline constexpr auto
RingmapMatrixRow::original_end_index() const noexcept -> base_index_type {
  return this->end_index_;
}

inline constexpr auto
RingmapMatrixRow::window_begin_index() const noexcept -> base_index_type {
  return this->window_begin_index_;
}

inline constexpr auto
RingmapMatrixRow::window_end_index() const noexcept -> base_index_type {
  return this->window_end_index_;
}

inline constexpr void RingmapMatrixRow::unsafe_set_original_begin_index(
    base_index_type begin_index) noexcept {
  this->begin_index_ = begin_index;
}

inline constexpr void RingmapMatrixRow::unsafe_set_original_end_index(
    base_index_type end_index) noexcept {
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

inline constexpr void
RingmapMatrixRow::copy_window_begin_end_indices(auto const &window) noexcept {
  using window_t = std::remove_cvref_t<decltype(window)>;

  if constexpr (std::is_same_v<window_t, RingmapMatrixWindowIndices> ||
                std::is_same_v<window_t, results::Window>) {
    this->window_begin_index_ = window.begin_index;
    this->window_end_index_ = window.end_index;
  } else if constexpr (std::is_same_v<window_t, RingmapMatrixRow>) {
    this->begin_index_ = window.begin_index_;
    this->end_index_ = window.end_index_;
  } else if constexpr (std::is_same_v<
                           window_t, RingmapMatrixRowAccessor<RingmapMatrix>> ||
                       std::is_same_v<window_t, RingmapMatrixRowAccessor<
                                                    const RingmapMatrix>>) {
    this->begin_index_ = window.original_begin_index();
    this->end_index_ = window.original_end_index();
  } else {
    static_assert(false, "copy_window_begin_end_indices only accepts a "
                         "RingmapMatrixWindowIndices, a results::Window, a "
                         "RingmapMatrixRow or a RingmapMatrixRowAccessor");
  }
}

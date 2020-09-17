#pragma once

#include "ringmap_matrix_traits.hpp"

#include <limits>
#include <vector>

struct RingmapMatrixRow : std::vector<ringmap_matrix::base_index_type> {
  using base_index_type = ringmap_matrix::base_index_type;
  using base_type = std::vector<base_index_type>;

  using base_type::base_type;
  using base_type::operator=;

  RingmapMatrixRow(const RingmapMatrixRow&) = default;
  RingmapMatrixRow(RingmapMatrixRow&&) = default;

  RingmapMatrixRow& operator=(const RingmapMatrixRow&) = default;
  RingmapMatrixRow& operator=(RingmapMatrixRow&&) = default;

  base_index_type begin_index = std::numeric_limits<base_index_type>::max();
  base_index_type end_index = std::numeric_limits<base_index_type>::min();
};

inline bool
operator<(RingmapMatrixRow const& lhs, RingmapMatrixRow const& rhs) noexcept {
  using base_type = RingmapMatrixRow::base_type;

  if (lhs.begin_index < rhs.begin_index) {
    return true;
  } else if (lhs.begin_index == rhs.begin_index) {
    if (lhs.end_index < rhs.end_index) {
      return true;
    } else if (lhs.end_index == rhs.end_index) {
      return static_cast<base_type const&>(lhs) <
             static_cast<base_type const&>(rhs);
    }
  }

  return false;
}

inline bool
operator==(RingmapMatrixRow const& lhs, RingmapMatrixRow const& rhs) noexcept {
  using base_type = RingmapMatrixRow::base_type;

  return lhs.begin_index == rhs.begin_index and
         lhs.end_index == rhs.end_index and
         static_cast<base_type const&>(lhs) ==
             static_cast<base_type const&>(rhs);
}

inline bool
operator!=(RingmapMatrixRow const& lhs, RingmapMatrixRow const& rhs) noexcept {
  return not(lhs == rhs);
}

inline bool
operator<=(RingmapMatrixRow const& lhs, RingmapMatrixRow const& rhs) noexcept {
  return lhs < rhs and lhs == rhs;
}

inline bool
operator>(RingmapMatrixRow const& lhs, RingmapMatrixRow const& rhs) noexcept {
  return not(lhs < rhs) and not(lhs == rhs);
}

inline bool
operator>=(RingmapMatrixRow const& lhs, RingmapMatrixRow const& rhs) noexcept {
  return not(lhs < rhs) and lhs == rhs;
}

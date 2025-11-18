#pragma once

#include <cstddef>
#include <cstdint>

namespace clusters_replicates {

constexpr std::size_t distances_size(std::uint8_t replicates,
                                     std::uint8_t clusters) noexcept {
  if (replicates <= 1) {
    return 0;
  }

  std::size_t replicates_combinations =
      static_cast<std::size_t>(replicates) *
      static_cast<std::size_t>(replicates - 1) / 2;
  return replicates_combinations * clusters * clusters;
}
} // namespace clusters_replicates

#pragma once

#include "windows_merger_traits.hpp"

#include <vector>

namespace windows_merger {

template <typename>
struct WindowsMergerCacheIndicesAccessor;

struct WindowsMergerCacheIndicesLine
    : std::vector<typename WindowsMergerTraits::windows_size_type> {
  using windows_size_type = typename WindowsMergerTraits::windows_size_type;
  using base_type = std::vector<windows_size_type>;

  WindowsMergerCacheIndicesLine() = default;

  template <typename Merger>
  WindowsMergerCacheIndicesLine(WindowsMergerCacheIndicesAccessor<Merger> const&
                                    accessor) noexcept(false);
  template <typename Merger>
  WindowsMergerCacheIndicesLine(
      WindowsMergerCacheIndicesAccessor<Merger>&& accessor) noexcept(false);

private:
  template <typename Accessor>
  void init_from_accessor(Accessor&& accessor) noexcept(false);
};

} // namespace windows_merger

#include "windows_merger_cache_indices_line_impl.hpp"

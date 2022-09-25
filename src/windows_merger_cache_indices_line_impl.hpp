#pragma once

#include "windows_merger_cache_indices_accessor.hpp"
#include "windows_merger_cache_indices_line.hpp"

#include <range/v3/algorithm.hpp>

namespace windows_merger {

template <typename Merger>
inline WindowsMergerCacheIndicesLine::WindowsMergerCacheIndicesLine(
    WindowsMergerCacheIndicesAccessor<Merger> const &accessor) noexcept(false) {
  init_from_accessor(accessor);
}

template <typename Merger>
inline WindowsMergerCacheIndicesLine::WindowsMergerCacheIndicesLine(
    WindowsMergerCacheIndicesAccessor<Merger> &&accessor) noexcept(false) {
  init_from_accessor(std::move(accessor));
}

template <typename Accessor>
inline void WindowsMergerCacheIndicesLine::init_from_accessor(
    Accessor &&accessor) noexcept(false) {
  base_type::resize(accessor.size());
  ranges::copy(accessor, ranges::begin(*this));
}

} // namespace windows_merger

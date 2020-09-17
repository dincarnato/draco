#pragma once

#include "windows_merger_traits.hpp"

#include <cassert>
#include <limits>
#include <type_traits>

namespace windows_merger {

template <typename Weight>
struct WindowsMergerWindowBaseWeightsAccessor {
  using weight_type = Weight;
  using decayed_weight_type = std::decay_t<weight_type>;
  using pointer = weight_type*;
  using reference = weight_type&;
  using iterator = pointer;
  using clusters_size_type = typename WindowsMergerTraits::clusters_size_type;

  WindowsMergerWindowBaseWeightsAccessor(pointer first, pointer last) noexcept
      : _begin(first), _end(last) {
    assert(first <= last);
  }

  constexpr iterator
  begin() const noexcept {
    return _begin;
  }

  constexpr iterator
  end() const noexcept {
    return _end;
  }

  constexpr clusters_size_type
  size() const noexcept {
    const auto size = _end - _begin;
    assert(size <= std::numeric_limits<clusters_size_type>::max());
    return static_cast<clusters_size_type>(size);
  }

  reference operator[](clusters_size_type index) const noexcept {
    assert(index < _end - _begin);
    return _begin[index];
  }

private:
  pointer _begin;
  pointer _end;
};

} // namespace windows_merger

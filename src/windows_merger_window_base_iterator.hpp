#pragma once

#include "windows_merger_traits.hpp"

#include <cassert>
#include <iterator>
#include <type_traits>

namespace windows_merger {

template <typename Merger> struct WindowsMergerWindowBaseIterator {
  using merger_type = std::decay_t<Merger>;
  using merger_pointer_type = std::remove_reference_t<Merger> *;
  using traits_type = WindowsMergerTraits;
  using windows_size_type = typename traits_type::windows_size_type;
  using bases_size_type = typename traits_type::bases_size_type;
  using signed_bases_size_type = std::conditional_t<
      std::is_same_v<bases_size_type, std::uint8_t>, std::int16_t,
      std::conditional_t<
          std::is_same_v<bases_size_type, std::uint16_t>, std::int32_t,
          std::conditional_t<std::is_same_v<bases_size_type, std::uint32_t>,
                             std::int64_t, std::ptrdiff_t>>>;
  using base_accessor = WindowsMergerWindowBaseAccessor<Merger>;
  using self = WindowsMergerWindowBaseIterator;

  using value_type = WindowsMergerWindowBase;
  using difference_type = signed_bases_size_type;
  using reference = base_accessor;
  using iterator_category = std::random_access_iterator_tag;

  WindowsMergerWindowBaseIterator() = default;
  explicit WindowsMergerWindowBaseIterator(
      Merger &merger, windows_size_type window_index,
      signed_bases_size_type base_index = signed_bases_size_type(0)) noexcept;

  reference operator*() const noexcept;
  reference operator[](difference_type index) const noexcept;

  self &operator++() noexcept;
  self operator++(int) noexcept;
  self &operator--() noexcept;
  self operator--(int) noexcept;

  self &operator+=(difference_type offset) noexcept;
  self &operator-=(difference_type offset) noexcept;

  // In-situ friend declarations
  friend self operator+(const self &lhs, difference_type offset) noexcept {
    return self(*lhs.merger, lhs.window_index,
                static_cast<bases_size_type>(lhs.base_index + offset));
  }

  friend self operator+(typename self::difference_type offset,
                        const self &lhs) noexcept {
    return lhs + offset;
  }

  friend self operator-(const self &lhs,
                        typename self::difference_type offset) noexcept {
    return self(*lhs.merger, lhs.window_index,
                static_cast<bases_size_type>(lhs.base_index - offset));
  }

  friend difference_type operator-(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    assert(lhs.window_index == rhs.window_index);

    return static_cast<difference_type>(lhs.base_index) -
           static_cast<difference_type>(rhs.base_index);
  }

  friend bool operator==(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    assert(lhs.window_index == rhs.window_index);
    return lhs.base_index == rhs.base_index;
  }

  friend bool operator!=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    assert(lhs.window_index == rhs.window_index);
    return lhs.base_index != rhs.base_index;
  }

  friend bool operator<(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    assert(lhs.window_index == rhs.window_index);
    return lhs.base_index < rhs.base_index;
  }

  friend bool operator<=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    assert(lhs.window_index == rhs.window_index);
    return lhs.base_index <= rhs.base_index;
  }

  friend bool operator>=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    assert(lhs.window_index == rhs.window_index);
    return lhs.base_index >= rhs.base_index;
  }

  friend bool operator>(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    assert(lhs.window_index == rhs.window_index);
    return lhs.base_index > rhs.base_index;
  }

private:
  merger_pointer_type merger;
  windows_size_type window_index;
  signed_bases_size_type base_index;
};

} // namespace windows_merger

#include "windows_merger_window_base_iterator_impl.hpp"

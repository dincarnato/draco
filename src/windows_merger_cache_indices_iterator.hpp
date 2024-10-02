#pragma once

#include "windows_merger_traits.hpp"

#include "nostd/type_traits.hpp"

namespace windows_merger {

template <typename> struct WindowsMergerCacheIndicesAccessor;

struct WindowsMergerCacheIndicesLine;

template <typename Merger> struct WindowsMergerCacheIndicesIterator {
  using merger_type = Merger;
  using decayed_merger_type = std::decay_t<Merger>;
  using merger_pointer_type = std::remove_reference_t<Merger> *;
  using traits_type = WindowsMergerTraits;
  using windows_size_type = typename traits_type::windows_size_type;
  using signed_windows_size_type = std::conditional_t<
      std::is_same_v<windows_size_type, std::uint8_t>, std::int16_t,
      std::conditional_t<
          std::is_same_v<windows_size_type, std::uint16_t>, std::int32_t,
          std::conditional_t<std::is_same_v<windows_size_type, std::uint32_t>,
                             std::int64_t, std::ptrdiff_t>>>;
  using accessor = WindowsMergerCacheIndicesAccessor<Merger>;
  using self = WindowsMergerCacheIndicesIterator;

  using value_type = WindowsMergerCacheIndicesLine;
  using difference_type = signed_windows_size_type;
  using reference = accessor;
  using iterator_category = std::random_access_iterator_tag;

  WindowsMergerCacheIndicesIterator() = default;
  explicit WindowsMergerCacheIndicesIterator(
      Merger &merger, signed_windows_size_type line_index =
                          signed_windows_size_type(0)) noexcept;
  explicit WindowsMergerCacheIndicesIterator(
      Merger &merger, windows_size_type line_index) noexcept;

  WindowsMergerCacheIndicesIterator(const self &) = default;
  WindowsMergerCacheIndicesIterator(self &&) = default;
  self &operator=(const self &) = default;
  self &operator=(self &&) = default;

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
    return self(*lhs.merger,
                static_cast<windows_size_type>(lhs.line_index + offset));
  }

  friend self operator+(typename self::difference_type offset,
                        const self &lhs) noexcept {
    return lhs + offset;
  }

  friend self operator-(const self &lhs,
                        typename self::difference_type offset) noexcept {
    return self(*lhs.merger,
                static_cast<windows_size_type>(lhs.line_index - offset));
  }

  friend difference_type operator-(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);

    return static_cast<difference_type>(lhs.line_index) -
           static_cast<difference_type>(rhs.line_index);
  }

  friend bool operator==(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.line_index == rhs.line_index;
  }

  friend bool operator!=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.line_index != rhs.line_index;
  }

  friend bool operator<(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.line_index < rhs.line_index;
  }

  friend bool operator<=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.line_index <= rhs.line_index;
  }

  friend bool operator>=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.line_index >= rhs.line_index;
  }

  friend bool operator>(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.line_index > rhs.line_index;
  }

private:
  merger_pointer_type merger;
  signed_windows_size_type line_index;
};

} // namespace windows_merger

#include "windows_merger_cache_indices_iterator_impl.hpp"

#pragma once

#include "windows_merger_traits.hpp"
#include "windows_merger_window.hpp"

#include <cassert>
#include <iterator>
#include <type_traits>

namespace windows_merger {

template <typename> struct WindowsMergerWindowAccessor;

template <typename Merger> struct WindowsMergerWindowsIterator {
  using merger_pointer = std::remove_reference_t<Merger> *;
  using traits_type = WindowsMergerTraits;
  using windows_size_type = typename traits_type::windows_size_type;
  using signed_windows_size_type = std::conditional_t<
      std::is_same_v<windows_size_type, std::uint8_t>, std::int16_t,
      std::conditional_t<
          std::is_same_v<windows_size_type, std::uint16_t>, std::int32_t,
          std::conditional_t<std::is_same_v<windows_size_type, std::uint32_t>,
                             std::int64_t, std::ptrdiff_t>>>;
  using accessor = WindowsMergerWindowAccessor<Merger>;
  using self = WindowsMergerWindowsIterator;

  using value_type = WindowsMergerWindow;
  using difference_type = signed_windows_size_type;
  using reference = accessor;
  using iterator_category = std::random_access_iterator_tag;

  WindowsMergerWindowsIterator() = default;
  explicit WindowsMergerWindowsIterator(
      Merger &merger, signed_windows_size_type window_index =
                          signed_windows_size_type(0)) noexcept;
  explicit WindowsMergerWindowsIterator(
      Merger &merger, windows_size_type window_index) noexcept;

  WindowsMergerWindowsIterator(const self &) = default;
  WindowsMergerWindowsIterator(self &&) = default;
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
                static_cast<windows_size_type>(lhs.window_index + offset));
  }

  friend self operator+(typename self::difference_type offset,
                        const self &lhs) noexcept {
    return lhs + offset;
  }

  friend self operator-(const self &lhs,
                        typename self::difference_type offset) noexcept {
    return self(*lhs.merger,
                static_cast<windows_size_type>(lhs.window_index - offset));
  }

  friend difference_type operator-(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);

    return static_cast<difference_type>(lhs.window_index) -
           static_cast<difference_type>(rhs.window_index);
  }

  friend bool operator==(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.window_index == rhs.window_index;
  }

  friend bool operator!=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.window_index != rhs.window_index;
  }

  friend bool operator<(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.window_index < rhs.window_index;
  }

  friend bool operator<=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.window_index <= rhs.window_index;
  }

  friend bool operator>=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.window_index >= rhs.window_index;
  }

  friend bool operator>(const self &lhs, const self &rhs) noexcept {
    assert(lhs.merger == rhs.merger);
    return lhs.window_index > rhs.window_index;
  }

private:
  merger_pointer merger;
  signed_windows_size_type window_index;
};

} // namespace windows_merger

#include "windows_merger_windows_iterator_impl.hpp"

#pragma once

#include "windows_merger_traits.hpp"

#include "nostd/type_traits.hpp"

#include <cassert>
#include <iterator>

namespace windows_merger {

template <typename> struct WindowsMergerWindowCoveragesAccessor;

template <typename Window> struct WindowsMergerWindowCoveragesIterator {
  using window_type = std::decay_t<Window>;
  using window_pointer_type = std::remove_reference_t<Window> *;
  using traits_type = WindowsMergerTraits;
  using bases_size_type = typename traits_type::bases_size_type;
  using coverage_type = typename traits_type::coverage_type;
  using signed_bases_size_type = std::conditional_t<
      std::is_same_v<bases_size_type, std::uint8_t>, std::int16_t,
      std::conditional_t<
          std::is_same_v<bases_size_type, std::uint16_t>, std::int32_t,
          std::conditional_t<std::is_same_v<bases_size_type, std::uint32_t>,
                             std::int64_t, std::ptrdiff_t>>>;
  using base_accessor = WindowsMergerWindowCoveragesAccessor<Window>;
  using self = WindowsMergerWindowCoveragesIterator;

  using value_type = coverage_type;
  using difference_type = signed_bases_size_type;
  using reference = std::conditional_t<
      std::is_rvalue_reference_v<Window>, coverage_type &&,
      std::conditional_t<std::is_const_v<Window>, coverage_type const &,
                         coverage_type &>>;
  using pointer =
      nostd::copy_const_t<std::remove_reference_t<Window>, coverage_type> *;
  using iterator_category = std::random_access_iterator_tag;

  WindowsMergerWindowCoveragesIterator() = default;
  explicit WindowsMergerWindowCoveragesIterator(
      Window &window,
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
    return self(*lhs.window,
                static_cast<bases_size_type>(lhs.base_index + offset));
  }

  friend self operator+(typename self::difference_type offset,
                        const self &lhs) noexcept {
    return lhs + offset;
  }

  friend self operator-(const self &lhs,
                        typename self::difference_type offset) noexcept {
    return self(*lhs.window,
                static_cast<bases_size_type>(lhs.base_index - offset));
  }

  friend difference_type operator-(const self &lhs, const self &rhs) noexcept {
    assert(lhs.window == rhs.window);

    return static_cast<difference_type>(lhs.base_index) -
           static_cast<difference_type>(rhs.base_index);
  }

  friend bool operator==(const self &lhs, const self &rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.base_index == rhs.base_index;
  }

  friend bool operator!=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.base_index != rhs.base_index;
  }

  friend bool operator<(const self &lhs, const self &rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.base_index < rhs.base_index;
  }

  friend bool operator<=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.base_index <= rhs.base_index;
  }

  friend bool operator>=(const self &lhs, const self &rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.base_index >= rhs.base_index;
  }

  friend bool operator>(const self &lhs, const self &rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.base_index > rhs.base_index;
  }

private:
  window_pointer_type window;
  signed_bases_size_type base_index;
};

} // namespace windows_merger

#include "windows_merger_window_coverages_iterator_impl.hpp"

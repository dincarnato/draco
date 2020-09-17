#pragma once

#include "windows_merger_traits.hpp"

#include "nostd/type_traits.hpp"

#include <cassert>
#include <iterator>
#include <range/v3/core.hpp>

namespace windows_merger {

template <typename>
struct WindowsMergerWindowWeightsAccessor;

template <typename Window>
struct WindowsMergerWindowWeightsIterator {
  using window_type = std::decay_t<Window>;
  using window_pointer_type = std::remove_reference_t<Window>*;
  using traits_type = WindowsMergerTraits;
  using bases_size_type = typename traits_type::bases_size_type;
  using weight_type = typename traits_type::weight_type;
  using clusters_size_type = typename traits_type::clusters_size_type;
  using signed_bases_size_type = std::conditional_t<
      std::is_same_v<bases_size_type, std::uint8_t>, std::int16_t,
      std::conditional_t<
          std::is_same_v<bases_size_type, std::uint16_t>, std::int32_t,
          std::conditional_t<std::is_same_v<bases_size_type, std::uint32_t>,
                             std::int64_t, std::ptrdiff_t>>>;
  using base_accessor = WindowsMergerWindowWeightsAccessor<Window>;
  using self = WindowsMergerWindowWeightsIterator;

  using value_type = weight_type;
  using difference_type = signed_bases_size_type;
  using reference = std::conditional_t<
      std::is_rvalue_reference_v<Window>, weight_type&&,
      std::conditional_t<std::is_const_v<Window>, weight_type const&,
                         weight_type&>>;
  using pointer =
      nostd::copy_const_t<std::remove_reference_t<Window>, weight_type>*;
  using iterator_category = ranges::random_access_iterator_tag;

  WindowsMergerWindowWeightsIterator() = default;
  explicit WindowsMergerWindowWeightsIterator(
      Window& window,
      signed_bases_size_type weight_index = signed_bases_size_type(0)) noexcept;

  reference operator*() const noexcept;
  reference operator[](difference_type index) const noexcept;

  self& operator++() noexcept;
  self operator++(int) noexcept;
  self& operator--() noexcept;
  self operator--(int) noexcept;

  self& operator+=(difference_type offset) noexcept;
  self& operator-=(difference_type offset) noexcept;

  // In-situ friend declarations
  friend self
  operator+(const self& lhs, difference_type offset) noexcept {
    return self(*lhs.window,
                static_cast<bases_size_type>(lhs.weight_index + offset));
  }

  friend self
  operator+(typename self::difference_type offset, const self& lhs) noexcept {
    return lhs + offset;
  }

  friend self
  operator-(const self& lhs, typename self::difference_type offset) noexcept {
    return self(*lhs.window,
                static_cast<bases_size_type>(lhs.weight_index - offset));
  }

  friend difference_type
  operator-(const self& lhs, const self& rhs) noexcept {
    assert(lhs.window == rhs.window);

    return static_cast<difference_type>(lhs.weight_index) -
           static_cast<difference_type>(rhs.weight_index);
  }

  friend bool
  operator==(const self& lhs, const self& rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.weight_index == rhs.weight_index;
  }

  friend bool
  operator!=(const self& lhs, const self& rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.weight_index != rhs.weight_index;
  }

  friend bool
  operator<(const self& lhs, const self& rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.weight_index < rhs.weight_index;
  }

  friend bool
  operator<=(const self& lhs, const self& rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.weight_index <= rhs.weight_index;
  }

  friend bool
  operator>=(const self& lhs, const self& rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.weight_index >= rhs.weight_index;
  }

  friend bool
  operator>(const self& lhs, const self& rhs) noexcept {
    assert(lhs.window == rhs.window);
    return lhs.weight_index > rhs.weight_index;
  }

private:
  window_pointer_type window;
  signed_bases_size_type weight_index;
  clusters_size_type n_clusters;
};

} // namespace windows_merger

#include "windows_merger_window_weights_iterator_impl.hpp"

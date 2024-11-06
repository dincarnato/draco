#pragma once

#include "windows_merger_traits.hpp"

#include "nostd/type_traits.hpp"

namespace windows_merger {

template <typename Window> struct WindowsMergerWindowCoveragesIterator;

template <typename Window> struct WindowsMergerWindowCoveragesAccessor {
  using window_type = Window;
  using decay_window_type = std::decay_t<Window>;
  using window_pointer_type = std::remove_reference_t<Window> *;
  using traits_type = WindowsMergerTraits;
  using bases_size_type = typename traits_type::bases_size_type;
  using coverage_type = typename traits_type::coverage_type;
  using reference = std::conditional_t<
      std::is_rvalue_reference_v<Window>, coverage_type &&,
      std::conditional_t<std::is_const_v<Window>, coverage_type const &,
                         coverage_type &>>;
  using iterator = WindowsMergerWindowCoveragesIterator<Window>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  WindowsMergerWindowCoveragesAccessor(Window &window) noexcept;

  bases_size_type size() const noexcept;

  reference front() const noexcept;
  reference back() const noexcept;
  reference operator[](bases_size_type index) const noexcept;

  iterator begin() const noexcept;
  iterator end() const noexcept;
  reverse_iterator rbegin() const noexcept;
  reverse_iterator rend() const noexcept;

private:
  window_pointer_type _window;
};

} // namespace windows_merger

#include "windows_merger_window_coverages_accessor_impl.hpp"

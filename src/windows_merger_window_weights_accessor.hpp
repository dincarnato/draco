#pragma once

#include "windows_merger_traits.hpp"

#include "nostd/type_traits.hpp"
#include <range/v3/core.hpp>

namespace windows_merger {

template <typename> struct WindowsMergerWindowWeightsIterator;

template <typename Window> struct WindowsMergerWindowWeightsAccessor {
  using window_type = Window;
  using decay_window_type = std::decay_t<Window>;
  using window_pointer_type = std::remove_reference_t<Window> *;
  using traits_type = WindowsMergerTraits;
  using bases_size_type = typename traits_type::bases_size_type;
  using weight_type = typename traits_type::weight_type;
  using reference = std::conditional_t<
      std::is_rvalue_reference_v<Window>, weight_type &&,
      std::conditional_t<std::is_const_v<Window>, weight_type const &,
                         weight_type &>>;
  using iterator = WindowsMergerWindowWeightsIterator<Window>;
  using reverse_iterator = ranges::reverse_iterator<iterator>;

  WindowsMergerWindowWeightsAccessor(Window &window) noexcept;

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

#include "windows_merger_window_weights_accessor_impl.hpp"

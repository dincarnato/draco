#pragma once

#include "windows_merger_traits.hpp"

#include "nostd/type_traits.hpp"

namespace windows_merger {

template <typename Merger> struct WindowsMergerWindowAccessorCoveragesAccessor {
  using merger_type = Merger;
  using merger_pointer_type = std::remove_reference_t<Merger> *;
  using traits_type = WindowsMergerTraits;
  using index_type = typename traits_type::windows_size_type;
  using bases_size_type = typename traits_type::bases_size_type;
  using coverage_type = typename traits_type::coverage_type;
  using reference = nostd::copy_const_t<Merger, coverage_type> &;
  using iterator = nostd::copy_const_t<Merger, coverage_type> *;

  WindowsMergerWindowAccessorCoveragesAccessor(
      Merger &merger, index_type window_index) noexcept;

  index_type index() const noexcept;
  bases_size_type size() const noexcept;

  reference front() const noexcept;
  reference back() const noexcept;
  reference operator[](bases_size_type index) const noexcept;

  iterator begin() const noexcept;
  iterator end() const noexcept;

private:
  merger_pointer_type merger;
  index_type window_index;
};

} // namespace windows_merger

#include "windows_merger_window_accessor_coverages_accessor_impl.hpp"

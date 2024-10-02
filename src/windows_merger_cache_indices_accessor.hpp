#pragma once

#include "windows_merger_traits.hpp"

#include "nostd/type_traits.hpp"

namespace windows_merger {

struct WindowsMergerCacheIndices;
struct WindowsMergerCacheIndicesLine;

template <typename Merger> struct WindowsMergerCacheIndicesAccessor {
  template <typename> friend struct WindowsMergerCacheIndicesAccessor;

  using merger_type = Merger;
  using decayed_merger_type = std::decay_t<Merger>;
  using merger_pointer_type = std::remove_reference_t<Merger> *;
  using traits_type = WindowsMergerTraits;
  using windows_size_type = typename traits_type::windows_size_type;
  using reference = nostd::copy_const_t<Merger, windows_size_type> &;
  using iterator = nostd::copy_const_t<Merger, windows_size_type> *;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using self = WindowsMergerCacheIndicesAccessor;

  WindowsMergerCacheIndicesAccessor(Merger &merger,
                                    windows_size_type line_index) noexcept;
  WindowsMergerCacheIndicesAccessor(WindowsMergerCacheIndicesAccessor const &) =
      default;
  WindowsMergerCacheIndicesAccessor(WindowsMergerCacheIndicesAccessor &&) =
      default;

  self const &operator=(
      WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices> const &rhs)
      const noexcept(false);
  self const &operator=(
      WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices const> const
          &rhs) const noexcept(false);
  self const &operator=(
      WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices &&> const
          &rhs) const noexcept(false);
  self const &operator=(
      WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices> &&rhs) const
      noexcept(false);
  self const &operator=(
      WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices const> &&rhs)
      const noexcept(false);
  self const &operator=(
      WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices &&> &&rhs)
      const noexcept(false);

  self const &operator=(WindowsMergerCacheIndicesLine const &line) const
      noexcept(false);
  self const &operator=(WindowsMergerCacheIndicesLine &&line) const
      noexcept(false);

  operator WindowsMergerCacheIndicesLine() const noexcept(false);

  windows_size_type index() const noexcept;
  windows_size_type size() const noexcept;
  bool empty() const noexcept;

  reference front() const noexcept;
  reference back() const noexcept;

  reference operator[](windows_size_type index) const noexcept;

  iterator begin() const noexcept;
  iterator end() const noexcept;
  reverse_iterator rbegin() const noexcept;
  reverse_iterator rend() const noexcept;

  windows_size_type &emplace_back(windows_size_type index) const
      noexcept(false);
  void push_back(windows_size_type index) const noexcept(false);

  void resize(windows_size_type value) const noexcept(false);
  void clear() const noexcept;

private:
  template <typename Accessor>
  void assign_from_accessor(Accessor &&accessor) const noexcept(false);

  template <typename Line>
  void assign_from_line(Line &&line) const noexcept(false);

  merger_pointer_type merger;
  windows_size_type line_index;
};

} // namespace windows_merger

#include "windows_merger_cache_indices_accessor_impl.hpp"

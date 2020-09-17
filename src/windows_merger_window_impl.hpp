#pragma once

#include "windows_merger_window.hpp"
#include "windows_merger_window_base_iterator.hpp"
#include "windows_merger_window_coverages_accessor.hpp"
#include "windows_merger_window_weights_accessor.hpp"
#include "windows_merger_windows_iterator.hpp"

#include "weighted_clusters.hpp"

#include <limits>
#include <range/v3/algorithm.hpp>
#include <range/v3/core.hpp>
#include <range/v3/iterator/insert_iterators.hpp>

namespace windows_merger {

inline WindowsMergerWindow::
    WindowsMergerWindow(clusters_size_type n_clusters) noexcept(
        std::is_nothrow_default_constructible_v<bases_type>)
    : n_clusters(n_clusters) {}

inline WindowsMergerWindow::WindowsMergerWindow(
    WindowsMergerWindowAccessor<WindowsMergerWindows> const&
        rhs) noexcept(false)
    : _bases(rhs.size()), _begin_index(rhs.begin_index()),
      n_clusters(rhs.clusters_size()) {
  ranges::copy(rhs, ranges::begin(_bases));
}

inline WindowsMergerWindow::WindowsMergerWindow(
    WindowsMergerWindowAccessor<WindowsMergerWindows>&& rhs) noexcept(false)
    : WindowsMergerWindow(rhs) {}

inline WindowsMergerWindow::WindowsMergerWindow(
    WindowsMergerWindowAccessor<WindowsMergerWindows const> const&
        rhs) noexcept(false)
    : _bases(rhs.size()), _begin_index(rhs.begin_index()),
      n_clusters(rhs.clusters_size()) {
  ranges::copy(rhs, ranges::begin(_bases));
}

inline WindowsMergerWindow::WindowsMergerWindow(
    WindowsMergerWindowAccessor<WindowsMergerWindows const>&&
        rhs) noexcept(false)
    : WindowsMergerWindow(rhs) {}

inline WindowsMergerWindow::WindowsMergerWindow(
    WindowsMergerWindowAccessor<WindowsMergerWindows&&> const&
        rhs) noexcept(false)
    : _bases(ranges::move_iterator(ranges::begin(rhs)),
             ranges::move_iterator(ranges::end(rhs))),
      _begin_index(rhs.begin_index()), n_clusters(rhs.clusters_size()) {}

inline WindowsMergerWindow::WindowsMergerWindow(
    WindowsMergerWindowAccessor<WindowsMergerWindows&&>&& rhs) noexcept(false)
    : WindowsMergerWindow(rhs) {}

inline auto
WindowsMergerWindow::
operator=(WindowsMergerWindowAccessor<WindowsMergerWindows> const&
              rhs) noexcept(false) -> self& {
  assign_from_accessor(rhs);
  return *this;
}

inline auto
WindowsMergerWindow::operator=(
    WindowsMergerWindowAccessor<WindowsMergerWindows>&& rhs) noexcept(false)
    -> self& {
  assign_from_accessor(std::move(rhs));
  return *this;
}

inline auto
WindowsMergerWindow::
operator=(WindowsMergerWindowAccessor<WindowsMergerWindows const> const&
              rhs) noexcept(false) -> self& {
  assign_from_accessor(rhs);
  return *this;
}

inline auto
WindowsMergerWindow::
operator=(WindowsMergerWindowAccessor<WindowsMergerWindows const>&&
              rhs) noexcept(false) -> self& {
  assign_from_accessor(std::move(rhs));
  return *this;
}

inline auto
WindowsMergerWindow::
operator=(WindowsMergerWindowAccessor<WindowsMergerWindows&&> const&
              rhs) noexcept(false) -> self& {
  assign_from_accessor(rhs);
  return *this;
}

inline auto
WindowsMergerWindow::operator=(
    WindowsMergerWindowAccessor<WindowsMergerWindows&&>&& rhs) noexcept(false)
    -> self& {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Accessor>
void
WindowsMergerWindow::assign_from_accessor(Accessor&& accessor) noexcept(false) {
  n_clusters = accessor.clusters_size();
  _begin_index = accessor.begin_index();
  _bases.clear();
  _bases.reserve(accessor.size());

  using merger_type = typename std::decay_t<Accessor>::merger_type;

  if constexpr (std::is_rvalue_reference_v<merger_type>)
    ranges::move(ranges::begin(accessor), ranges::end(accessor),
                 ranges::back_inserter(_bases));
  else
    ranges::copy(ranges::begin(accessor), ranges::end(accessor),
                 ranges::back_inserter(_bases));
}

inline auto
WindowsMergerWindow::begin_index() const noexcept -> bases_size_type {
  return _begin_index;
}

inline auto
WindowsMergerWindow::end_index() const noexcept -> bases_size_type {
  assert(_bases.size() <=
         static_cast<std::size_t>(std::numeric_limits<bases_size_type>::max() -
                                  _begin_index));
  return static_cast<bases_size_type>(_begin_index + _bases.size());
}

inline auto
WindowsMergerWindow::size() const noexcept -> bases_size_type {
  assert(_bases.size() <= std::numeric_limits<bases_size_type>::max());
  return static_cast<bases_size_type>(_bases.size());
}

inline bool
WindowsMergerWindow::empty() const noexcept {
  return _bases.empty();
}

inline auto
WindowsMergerWindow::clusters_size() const noexcept -> clusters_size_type {
  return n_clusters;
}

inline auto
WindowsMergerWindow::coverages() noexcept -> coverages_accessor {
  return coverages_accessor(*this);
}

inline auto
WindowsMergerWindow::coverages() const noexcept -> const_coverages_accessor {
  return const_coverages_accessor(*this);
}

inline auto
WindowsMergerWindow::linear_weights() noexcept -> linear_weights_accessor {
  return linear_weights_accessor(*this);
}

inline auto
WindowsMergerWindow::linear_weights() const noexcept
    -> const_linear_weights_accessor {
  return const_linear_weights_accessor(*this);
}

inline void
WindowsMergerWindow::set_begin_index(bases_size_type begin_index) noexcept {
  _begin_index = begin_index;
}

inline auto
WindowsMergerWindow::front() const noexcept -> base_type const& {
  assert(not _bases.empty());
  return _bases.front();
}

inline auto
WindowsMergerWindow::back() const noexcept -> base_type const& {
  assert(not _bases.empty());
  return _bases.back();
}

inline auto
WindowsMergerWindow::front() noexcept -> base_type& {
  assert(not _bases.empty());
  return _bases.front();
}

inline auto
WindowsMergerWindow::back() noexcept -> base_type& {
  assert(not _bases.empty());
  return _bases.back();
}

inline auto
    WindowsMergerWindow::begin() &
    noexcept -> iterator {
  return _bases.begin();
}

inline auto
    WindowsMergerWindow::end() &
    noexcept -> iterator {
  return _bases.end();
}

inline auto
    WindowsMergerWindow::begin() const
    & noexcept -> const_iterator {
  return _bases.begin();
}

inline auto
    WindowsMergerWindow::end() const
    & noexcept -> const_iterator {
  return _bases.end();
}

inline auto
    WindowsMergerWindow::begin() &&
    noexcept -> move_iterator {
  return ranges::move_iterator(_bases.begin());
}

inline auto
    WindowsMergerWindow::end() &&
    noexcept -> move_iterator {
  return ranges::move_iterator(_bases.end());
}

inline auto
    WindowsMergerWindow::rbegin() &
    noexcept -> reverse_iterator {
  return _bases.rbegin();
}

inline auto
    WindowsMergerWindow::rend() &
    noexcept -> reverse_iterator {
  return _bases.rend();
}

inline auto
    WindowsMergerWindow::rbegin() const
    & noexcept -> const_reverse_iterator {
  return _bases.rbegin();
}

inline auto
    WindowsMergerWindow::rend() const
    & noexcept -> const_reverse_iterator {
  return _bases.rend();
}

inline auto
    WindowsMergerWindow::rbegin() &&
    noexcept -> move_reverse_iterator {
  return ranges::move_iterator(_bases.rbegin());
}

inline auto
    WindowsMergerWindow::rend() &&
    noexcept -> move_reverse_iterator {
  return ranges::move_iterator(_bases.rend());
}

inline auto WindowsMergerWindow::operator[](bases_size_type index) const
    noexcept -> base_type const& {
  assert(index < _bases.size());
  return _bases[index];
}

inline auto WindowsMergerWindow::operator[](bases_size_type index) noexcept
    -> base_type& {
  assert(index < _bases.size());
  return _bases[index];
}

template <typename... Args>
inline auto
WindowsMergerWindow::emplace_back(Args&&... args) noexcept(false)
    -> base_type& {
  auto&& new_base = _bases.emplace_back(std::forward<Args>(args)...);
  assert(new_base.clusters_size() == n_clusters);
  return new_base;
}

template <typename Base>
inline void
WindowsMergerWindow::push_back(Base&& base) noexcept(false) {
  assert(base.clusters_size() == n_clusters);
  _bases.push_back(std::forward<Base>(base));
}

} // namespace windows_merger

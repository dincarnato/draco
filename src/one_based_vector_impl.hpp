#pragma once

#include "one_based_vector.hpp"

#include <algorithm>
#include <functional>
#include <limits>

template <typename T, typename Allocator>
one_based_vector<T, Allocator>::one_based_vector(
    const Allocator& allocator) noexcept
    : vec(allocator) {}

template <typename T, typename Allocator>
one_based_vector<T, Allocator>::one_based_vector(size_type count,
                                                 const T& value,
                                                 const Allocator& allocator)
    : vec([&] {
        if (count == std::numeric_limits<size_type>::max())
          throw std::overflow_error("size too big");

        zero_based_type vec(count + 1, allocator);
        std::fill(std::next(std::begin(vec)), std::end(vec), value);
        return vec;
      }()) {}

template <typename T, typename Allocator>
one_based_vector<T, Allocator>::one_based_vector(size_type count,
                                                 const Allocator& allocator)
    : vec([&] {
        if (count == std::numeric_limits<size_type>::max())
          throw std::overflow_error("size too big");

        return zero_based_type(count + 1, allocator);
      }()) {}

template <typename T, typename Allocator>
template <typename InputIt>
one_based_vector<T, Allocator>::one_based_vector(InputIt first, InputIt last,
                                                 const Allocator& alloc)
    : vec([&] {
        auto size = std::distance(first, last);
        if (size == std::numeric_limits<size_type>::max())
          throw std::overflow_error("size too big");
        else if (size <= 0)
          return zero_based_type(alloc);
        else {
          zero_based_type vec(size + 1, alloc);
          std::copy(std::move(first), std::move(last),
                    std::next(std::begin(vec)));

          return vec;
        }
      }()) {}

template <typename T, typename Allocator>
one_based_vector<T, Allocator>::one_based_vector(const one_based_vector& other,
                                                 const Allocator& alloc)
    : vec(other, alloc) {}

template <typename T, typename Allocator>
one_based_vector<T, Allocator>::one_based_vector(one_based_vector&& other,
                                                 const Allocator& alloc)
    : vec(std::move(other), alloc) {}

template <typename T, typename Allocator>
one_based_vector<T, Allocator>::one_based_vector(std::initializer_list<T> init,
                                                 const Allocator& alloc)
    : vec([&] {
        auto size = init.size();
        if (size == std::numeric_limits<size_type>::max())
          throw std::overflow_error("size too big");
        else if (size == 0)
          return zero_based_type(alloc);
        else {
          zero_based_type vec(size + 1, alloc);
          std::move(std::begin(init), std::end(init),
                    std::next(std::begin(vec)));

          return vec;
        }
      }()) {}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::operator=(std::initializer_list<T> ilist)
    -> one_based_vector& {
  auto size = ilist.size();
  if (size == std::numeric_limits<size_type>::max())
    throw std::overflow_error("size too big");

  vec.resize(size + 1);
  std::move(std::begin(ilist), std::end(ilist), std::next(std::begin(vec)));
  return *this;
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::assign(size_type count,
                                       const value_type& value) {
  vec.assign(count, value);
}

template <typename T, typename Allocator>
template <typename InputIt>
void
one_based_vector<T, Allocator>::assign(InputIt first, InputIt last) {
  auto size = std::distance(first, last);

  if (size == std::numeric_limits<size_type>::max())
    throw std::overflow_error("interval too big");
  else if (size <= 0)
    vec.clear();
  else {
    vec.resize(size + 1);
    std::copy(first, last, std::next(std::begin(vec)));
  }
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::assign(
    std::initializer_list<value_type> ilist) {
  auto size = ilist.size();
  if (size == std::numeric_limits<size_type>::max())
    throw std::overflow_error("too many elements");

  vec.resize(size + 1);
  std::move(std::begin(ilist), std::end(ilist), std::next(std::begin(vec)));
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::get_allocator() const -> allocator_type {
  return vec.get_allocator();
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::at(size_type pos) -> reference {
  if (pos == std::numeric_limits<size_type>::max())
    throw std::overflow_error("value too high");

  return vec.at(pos + 1);
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::at(size_type pos) const -> const_reference {
  if (pos == std::numeric_limits<size_type>::max())
    throw std::overflow_error("value too high");

  return vec.at(pos + 1);
}

template <typename T, typename Allocator>
auto one_based_vector<T, Allocator>::operator[](size_type pos) -> reference {
  if (pos == std::numeric_limits<size_type>::max())
    throw std::overflow_error("value too high");

  return vec[pos + 1];
}

template <typename T, typename Allocator>
auto one_based_vector<T, Allocator>::operator[](size_type pos) const
    -> const_reference {
  if (pos == std::numeric_limits<size_type>::max())
    throw std::overflow_error("value too high");

  return vec[pos + 1];
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::front() -> reference {
  return vec[1];
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::front() const -> const_reference {
  return vec[1];
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::back() -> reference {
  return vec.back();
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::back() const -> const_reference {
  return vec.back();
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::data() noexcept -> value_type* {
  return vec.data() + 1;
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::data() const noexcept -> const value_type* {
  return vec.data() + 1;
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::begin() noexcept -> iterator {
  if (vec.size() < 2)
    return vec.end();
  else
    return std::next(vec.begin());
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::begin() const noexcept -> const_iterator {
  if (vec.size() < 2)
    return vec.end();
  else
    return std::next(vec.begin());
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::cbegin() const noexcept -> const_iterator {
  if (vec.size() < 2)
    return vec.cend();
  else
    return std::next(vec.cbegin());
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::end() noexcept -> iterator {
  return vec.end();
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::end() const noexcept -> const_iterator {
  return vec.end();
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::cend() const noexcept -> const_iterator {
  return vec.cend();
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::rbegin() noexcept -> reverse_iterator {
  return vec.rbegin();
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::rbegin() const noexcept
    -> const_reverse_iterator {
  return vec.rbegin();
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::crbegin() const noexcept
    -> const_reverse_iterator {
  return vec.crbegin();
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::rend() noexcept -> reverse_iterator {
  if (vec.size() < 2)
    return vec.rbegin();
  else
    return std::prev(vec.rend());
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::rend() const noexcept
    -> const_reverse_iterator {
  if (vec.size() < 2)
    return vec.rbegin();
  else
    return std::prev(vec.rend());
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::crend() const noexcept
    -> const_reverse_iterator {
  if (vec.size() < 2)
    return vec.crbegin();
  else
    return std::prev(vec.crend());
}

template <typename T, typename Allocator>
[[nodiscard]] bool
one_based_vector<T, Allocator>::empty() const noexcept {
  return vec.size() < 2;
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::size() const noexcept -> size_type {
  auto size = vec.size();
  if (size < 2)
    return 0;
  else
    return size - 1;
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::max_size() const noexcept -> size_type {
  auto size = vec.max_size();
  if (size < 2)
    return 0;
  else
    return size - 1;
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::reserve(size_type new_cap) {
  if (new_cap == std::numeric_limits<size_type>::max())
    throw std::overflow_error("value too high");
  else if (new_cap > 0)
    vec.reserve(new_cap + 1);
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::capacity() const noexcept -> size_type {
  auto size = vec.capacity();
  if (size < 2)
    return 0;
  else
    return size - 1;
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::shrink_to_fit() {
  vec.shrink_to_fit();
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::clear() noexcept {
  vec.clear();
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::insert(const_iterator pos, const T& value)
    -> iterator {
  return vec.insert(std::move(pos), value);
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::insert(const_iterator pos, T&& value)
    -> iterator {
  return vec.insert(std::move(pos), std::move(value));
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::insert(const_iterator pos, size_type count,
                                       const T& value) -> iterator {
  return vec.insert(std::move(pos), std::move(count), value);
}

template <typename T, typename Allocator>
template <typename InputIt>
auto
one_based_vector<T, Allocator>::insert(const_iterator pos, InputIt first,
                                       InputIt last) -> iterator {
  return vec.insert(std::move(pos), std::move(first), std::move(last));
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::insert(const_iterator pos,
                                       std::initializer_list<T> ilist)
    -> iterator {
  return vec.insert(std::move(pos), std::move(ilist));
}

template <typename T, typename Allocator>
template <typename... Args>
auto
one_based_vector<T, Allocator>::emplace(const_iterator pos, Args&&... args)
    -> iterator {
  return vec.emplace(std::move(pos), std::forward<Args>(args)...);
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::erase(const_iterator pos) -> iterator {
  return vec.erase(std::move(pos));
}

template <typename T, typename Allocator>
auto
one_based_vector<T, Allocator>::erase(const_iterator first, const_iterator last)
    -> iterator {
  return vec.erase(std::move(first), std::move(last));
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::push_back(const T& value) {
  return vec.push_back(value);
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::push_back(T&& value) {
  return vec.push_back(std::move(value));
}

template <typename T, typename Allocator>
template <typename... Args>
auto
one_based_vector<T, Allocator>::emplace_back(Args&&... args) -> reference {
  return vec.emplace_back(std::forward<Args>(args)...);
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::pop_back() {
  vec.pop_back();
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::resize(size_type count) {
  if (count == std::numeric_limits<size_type>::max())
    throw std::overflow_error("value too high");
  else if (count == 0)
    vec.resize(0);
  else
    vec.resize(count + 1);
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::resize(size_type count,
                                       const value_type& value) {
  if (count == std::numeric_limits<size_type>::max())
    throw std::overflow_error("value too high");
  else if (count == 0)
    vec.resize(0);
  else {
    auto cur_size = std::min(vec.size(), static_cast<size_type>(1));
    vec.resize(count + 1);
    std::fill(std::next(std::begin(vec), cur_size), std::end(vec), value);
  }
}

template <typename T, typename Allocator>
void
one_based_vector<T, Allocator>::swap(one_based_vector& other) noexcept(
    noexcept(std::declval<zero_based_type&>().swap(other.vec))) {
  vec.swap(other.vec);
}

template <typename T, class Allocator>
bool
one_based_vector<T, Allocator>::operator==(const one_based_vector& rhs) {
  return std::equal(cbegin(), cend(), rhs.cbegin(), rhs.cend());
}

template <typename T, class Allocator>
bool
one_based_vector<T, Allocator>::operator!=(const one_based_vector& rhs) {
  return std::equal(cbegin(), cend(), rhs.cbegin(), rhs.cend(),
                    std::not_equal_to<T>());
}

template <typename T, class Allocator>
bool
one_based_vector<T, Allocator>::operator<(const one_based_vector& rhs) {
  return std::equal(cbegin(), cend(), rhs.cbegin(), rhs.cend(), std::less<T>());
}

template <typename T, class Allocator>
bool
one_based_vector<T, Allocator>::operator<=(const one_based_vector& rhs) {
  return std::equal(cbegin(), cend(), rhs.cbegin(), rhs.cend(),
                    std::less_equal<T>());
}

template <typename T, class Allocator>
bool
one_based_vector<T, Allocator>::operator>(const one_based_vector& rhs) {
  return std::equal(cbegin(), cend(), rhs.cbegin(), rhs.cend(),
                    std::greater<T>());
}

template <typename T, class Allocator>
bool
one_based_vector<T, Allocator>::operator>=(const one_based_vector& rhs) {
  return std::equal(cbegin(), cend(), rhs.cbegin(), rhs.cend(),
                    std::greater_equal<T>());
}

template <typename T, class Allocator>
auto
one_based_vector<T, Allocator>::zero_based() noexcept -> zero_based_type& {
  return vec;
}

template <typename T, class Allocator>
auto
one_based_vector<T, Allocator>::zero_based() const noexcept
    -> const zero_based_type& {
  return vec;
}

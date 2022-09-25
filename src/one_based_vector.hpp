#pragma once

#include <memory>
#include <vector>

template <typename T, typename Allocator = std::allocator<T>>
struct one_based_vector {
  using value_type = T;
  using allocator_type = Allocator;
  using zero_based_type = std::vector<T, Allocator>;
  using size_type = typename zero_based_type::size_type;
  using difference_type = typename zero_based_type::difference_type;
  using reference = typename zero_based_type::reference;
  using const_reference = typename zero_based_type::const_reference;
  using pointer = typename zero_based_type::pointer;
  using const_pointer = typename zero_based_type::const_pointer;
  using iterator = typename zero_based_type::iterator;
  using const_iterator = typename zero_based_type::const_iterator;
  using reverse_iterator = typename zero_based_type::reverse_iterator;
  using const_reverse_iterator =
      typename zero_based_type::const_reverse_iterator;

  one_based_vector() noexcept(noexcept(Allocator())) = default;
  explicit one_based_vector(const Allocator &alloc) noexcept;
  explicit one_based_vector(size_type count, const T &value,
                            const Allocator &alloc = Allocator());
  explicit one_based_vector(size_type count,
                            const Allocator &alloc = Allocator());
  template <typename InputIt>
  one_based_vector(InputIt first, InputIt last,
                   const Allocator &alloc = Allocator());
  one_based_vector(const one_based_vector &other) = default;
  one_based_vector(const one_based_vector &other, const Allocator &alloc);
  one_based_vector(one_based_vector &&other) noexcept = default;
  one_based_vector(one_based_vector &&other, const Allocator &alloc);
  one_based_vector(std::initializer_list<T> init,
                   const Allocator &alloc = Allocator());

  one_based_vector &operator=(const one_based_vector &other) = default;
  one_based_vector &operator=(one_based_vector &&other) noexcept(noexcept(
      std::declval<zero_based_type &>() = std::move(other.vec))) = default;
  one_based_vector &operator=(std::initializer_list<T> ilist);

  void assign(size_type count, const T &value);
  template <typename InputIt> void assign(InputIt first, InputIt last);
  void assign(std::initializer_list<T> ilist);

  Allocator get_allocator() const;

  reference at(size_type pos);
  const_reference at(size_type pos) const;

  reference operator[](size_type pos);
  const_reference operator[](size_type pos) const;

  reference front();
  const_reference front() const;

  reference back();
  const_reference back() const;

  T *data() noexcept;
  const T *data() const noexcept;

  iterator begin() noexcept;
  const_iterator begin() const noexcept;
  const_iterator cbegin() const noexcept;

  iterator end() noexcept;
  const_iterator end() const noexcept;
  const_iterator cend() const noexcept;

  reverse_iterator rbegin() noexcept;
  const_reverse_iterator rbegin() const noexcept;
  const_reverse_iterator crbegin() const noexcept;

  reverse_iterator rend() noexcept;
  const_reverse_iterator rend() const noexcept;
  const_reverse_iterator crend() const noexcept;

  [[nodiscard]] bool empty() const noexcept;
  size_type size() const noexcept;
  size_type max_size() const noexcept;
  void reserve(size_type new_cap);
  size_type capacity() const noexcept;
  void shrink_to_fit();
  void clear() noexcept;

  iterator insert(const_iterator pos, const T &value);
  iterator insert(const_iterator pos, T &&value);
  iterator insert(const_iterator pos, size_type count, const T &value);
  template <typename InputIt>
  iterator insert(const_iterator pos, InputIt first, InputIt last);
  iterator insert(const_iterator pos, std::initializer_list<T> ilist);

  template <typename... Args>
  iterator emplace(const_iterator pos, Args &&...args);

  iterator erase(const_iterator pos);
  iterator erase(const_iterator first, const_iterator last);

  void push_back(const T &value);
  void push_back(T &&value);

  template <typename... Args> reference emplace_back(Args &&...args);
  void pop_back();
  void resize(size_type count);
  void resize(size_type count, const value_type &value);
  void swap(one_based_vector &other) noexcept(
      noexcept(std::declval<zero_based_type &>().swap(other.vec)));

  bool operator==(const one_based_vector &rhs);
  bool operator!=(const one_based_vector &rhs);
  bool operator<(const one_based_vector &rhs);
  bool operator<=(const one_based_vector &rhs);
  bool operator>(const one_based_vector &rhs);
  bool operator>=(const one_based_vector &rhs);

  zero_based_type &zero_based() noexcept;
  const zero_based_type &zero_based() const noexcept;

private:
  zero_based_type vec;
};

template <typename InputIt,
          typename Alloc = std::allocator<
              typename std::iterator_traits<InputIt>::value_type>>
one_based_vector(InputIt, InputIt, Alloc = Alloc())
    -> one_based_vector<typename std::iterator_traits<InputIt>::value_type,
                        Alloc>;

#include "one_based_vector_impl.hpp"

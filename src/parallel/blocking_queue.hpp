#pragma once

#include "defaults.hpp"

#include <atomic>
#include <condition_variable>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>

namespace parallel {

template <typename T> struct blocking_queue_entry {
  template <typename, typename> friend class blocking_queue;

  using entry_type = T;

  template <typename... Args>
  blocking_queue_entry(Args &&...args) noexcept(
      noexcept(T(std::forward<Args>(args)...)));

  blocking_queue_entry(const blocking_queue_entry &) = delete;
  blocking_queue_entry(blocking_queue_entry &&) = delete;
  blocking_queue_entry &operator=(const blocking_queue_entry &) = delete;
  blocking_queue_entry &operator=(blocking_queue_entry &&) = delete;

#ifndef NDEBUG
  ~blocking_queue_entry() noexcept;
#endif

private:
  T entry;
  blocking_queue_entry *next_entry;
  blocking_queue_entry *prev_entry;
};

template <typename T, typename Alloc = std::allocator<blocking_queue_entry<T>>>
struct blocking_queue : private Alloc {
  using value_type = T;
  using reference = T &;
  using const_reference = const T &;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using allocator_type = Alloc;
  using entry_type = blocking_queue_entry<T>;

  blocking_queue() noexcept(noexcept(Alloc()));
  explicit blocking_queue(const allocator_type &alloc) noexcept;
  explicit blocking_queue(std::size_t max_size,
                          const allocator_type &alloc = Alloc()) noexcept;
  template <typename InputIterator>
  blocking_queue(InputIterator first, InputIterator last,
                 const allocator_type &alloc = Alloc());
  template <typename InputIterator>
  blocking_queue(InputIterator first, InputIterator last, std::size_t max_size,
                 const allocator_type &alloc = Alloc());
  blocking_queue(const blocking_queue &) = delete;
  blocking_queue(blocking_queue &&src) noexcept(
      std::is_nothrow_move_constructible_v<Alloc>);
  blocking_queue(blocking_queue &&src, const allocator_type &alloc) noexcept(
      std::is_nothrow_copy_constructible_v<Alloc>);
  blocking_queue &operator=(const blocking_queue &) = delete;
  blocking_queue &operator=(blocking_queue &&) noexcept(
      std::is_nothrow_move_assignable_v<Alloc>);

  ~blocking_queue() noexcept(std::is_nothrow_destructible_v<T>);

  template <typename U> void push(U &&value);

  template <typename... Args> void emplace(Args &&...args);

  std::optional<T> pop() noexcept(std::is_nothrow_move_constructible_v<T> and
                                  std::is_nothrow_destructible_v<T>);
  std::optional<T>
  try_pop() noexcept(std::is_nothrow_move_constructible_v<T> and
                     std::is_nothrow_destructible_v<T>);
  void finish() noexcept;
  bool finished() noexcept;

private:
  entry_type *head;
  std::atomic<entry_type *> tail;
  std::atomic<size_type> _size;
  std::size_t _max_size = std::numeric_limits<std::size_t>::max();
  std::atomic_bool pushing = false;
  std::atomic_bool popping = false;
  std::atomic_bool _finished = false;
  std::mutex mxEmpty, mxFull;
  std::condition_variable cvPop, cvPush;

  void moveDataFrom(blocking_queue &&other) noexcept;
  template <typename InputIterator>
  void initWithIterators(InputIterator first, InputIterator last);
};

template <typename InputIterator, typename Alloc>
blocking_queue(InputIterator, InputIterator, const Alloc &)
    -> blocking_queue<typename InputIterator::value_type, std::decay_t<Alloc>>;

} // namespace parallel

#include "blocking_queue_impl.hpp"

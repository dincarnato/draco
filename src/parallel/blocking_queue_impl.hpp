#pragma once

#include "blocking_queue.hpp"
#include <cassert>

namespace parallel {

template <typename T>
template <typename... Args>
blocking_queue_entry<T>::blocking_queue_entry(Args&&... args) noexcept(
    noexcept(T(std::forward<Args>(args)...)))
    : entry(std::forward<Args>(args)...), next_entry(nullptr) {}

#ifndef NDEBUG
template <typename T>
blocking_queue_entry<T>::~blocking_queue_entry() noexcept {
  /* We don't want to destroy an entry if it is not the last, or there is a high
   * risk of leaking. When not in debug, destroing an entry with a next_entry
   * different from nullptr is undefined behaviour */
  assert(next_entry == nullptr);
  assert(prev_entry == nullptr);
}
#endif

template <typename T, typename Alloc>
blocking_queue<T, Alloc>::blocking_queue() noexcept(noexcept(Alloc()))
    : Alloc(), head(nullptr), tail(nullptr), _size(0) {}

template <typename T, typename Alloc>
blocking_queue<T, Alloc>::blocking_queue(const allocator_type& alloc) noexcept
    : Alloc(alloc), head(nullptr), tail(nullptr), _size(0) {}

template <typename T, typename Alloc>
blocking_queue<T, Alloc>::blocking_queue(std::size_t max_size,
                                         const allocator_type& alloc) noexcept
    : Alloc(alloc), head(nullptr), tail(nullptr), _size(0),
      _max_size(max_size) {}

template <typename T, typename Alloc>
template <typename InputIterator>
void
blocking_queue<T, Alloc>::initWithIterators(InputIterator first,
                                            InputIterator last) {
  static_assert(std::is_same_v<typename InputIterator::value_type, T>,
                "InputIterator must have a value_type with the same type as T");

  std::size_t _size = 0;
  head = nullptr;
  entry_type* tail = nullptr;

  if (first != last) {
    entry_type* head = std::allocator_traits<Alloc>::allocate(*this, 1);
    std::allocator_traits<Alloc>::construct(*this, head, *first++);
    tail = head;
    ++_size;

    for (; first != last; ++first, ++_size) {
      entry_type* new_entry = std::allocator_traits<Alloc>::allocate(*this, 1);
      std::allocator_traits<Alloc>::construct(*this, new_entry, *first);

      head->prev_entry.store(new_entry, std::memory_order_release);
      new_entry->next_entry.store(std::exchange(head, new_entry),
                                  std::memory_order_release);
    }
  }

  this->tail.store(tail, std::memory_order_release);
  this->_size.store(_size, std::memory_order_release);
  head->prev_entry = nullptr;
}

template <typename T, typename Alloc>
template <typename InputIterator>
blocking_queue<T, Alloc>::blocking_queue(InputIterator first,
                                         InputIterator last,
                                         const allocator_type& a)
    : Alloc(a) {
  initWithIterators(std::move(first), std::move(last));
}

template <typename T, typename Alloc>
template <typename InputIterator>
blocking_queue<T, Alloc>::blocking_queue(InputIterator first,
                                         InputIterator last,
                                         std::size_t max_size,
                                         const allocator_type& a)
    : Alloc(a), _max_size(max_size) {
  initWithIterators(std::move(first), std::move(last));
}

template <typename T, typename Alloc>
blocking_queue<T, Alloc>::blocking_queue(blocking_queue&& other) noexcept(
    std::is_nothrow_move_constructible_v<Alloc>)
    : Alloc(static_cast<Alloc&&>(std::move(other))) {
  moveDataFrom(std::move(other));
}

template <typename T, typename Alloc>
blocking_queue<T, Alloc>::blocking_queue(
    blocking_queue&& other,
    const allocator_type&
        alloc) noexcept(std::is_nothrow_copy_constructible_v<Alloc>)
    : Alloc(alloc) {
  moveDataFrom(std::move(other));
}

template <typename T, typename Alloc>
void
blocking_queue<T, Alloc>::moveDataFrom(blocking_queue&& other) noexcept {
  for (bool pushing = false; not other.pushing.compare_exchange_weak(
           pushing, true, std::memory_order_acq_rel);)
    pushing = false;

  for (bool popping = false; not other.popping.compare_exchange_weak(
           popping, true, std::memory_order_acq_rel);)
    popping = false;

  head = std::exchange(other.head, nullptr);
  _size.store(other._size.exchange(static_cast<std::size_t>(0),
                                   std::memory_order_acq_rel),
              std::memory_order_release);
  tail.store(other.tail.exchange(nullptr, std::memory_order_acq_rel),
             std::memory_order_release);

  other.cvPop.notify_all();
  other.cvPush.notify_all();
}

template <typename T, typename Alloc>
auto
blocking_queue<T, Alloc>::operator=(blocking_queue&& other) noexcept(
    std::is_nothrow_move_assignable_v<Alloc>) -> blocking_queue& {
  static_cast<Alloc&>(*this) = static_cast<Alloc&&>(std::move(other));
  moveDataFrom(std::move(other));
  return *this;
}

template <typename T, typename Alloc>
blocking_queue<T, Alloc>::~blocking_queue() noexcept(
    std::is_nothrow_destructible_v<T>) {
  _finished.store(true, std::memory_order_release);

  for (entry_type* entry = std::exchange(head, nullptr); entry != nullptr;) {

#ifndef NDEBUG
    entry->prev_entry = nullptr;
    auto next_entry = std::exchange(entry->next_entry, nullptr);
#else
    auto next_entry = std::exchange(entry->next_entry, nullptr);
#endif
    std::allocator_traits<Alloc>::destroy(*this, entry);
    std::allocator_traits<Alloc>::deallocate(*this, entry, 1);

    entry = next_entry;
  }
  tail.store(nullptr, std::memory_order_release);
  cvPush.notify_all();
}

template <typename T, typename Alloc>
void
blocking_queue<T, Alloc>::finish() noexcept {
  while (not mxEmpty.try_lock())
    ;
  _finished.store(true, std::memory_order_release);
  cvPush.notify_all();
  mxEmpty.unlock();
}

template <typename T, typename Alloc>
bool
blocking_queue<T, Alloc>::finished() noexcept {
  return _finished.load(std::memory_order_acquire);
}

template <typename T, typename Alloc>
template <typename U>
void
blocking_queue<T, Alloc>::push(U&& value) {
  static_assert(std::is_same_v<std::decay_t<U>, T>,
                "U must decay to the same type of T");
  emplace(std::forward<U>(value));
}

template <typename T, typename Alloc>
template <typename... Args>
void
blocking_queue<T, Alloc>::emplace(Args&&... args) {
  static_assert(std::is_constructible_v<T, Args...>,
                "T must be constructible using Args...");

  entry_type* new_entry = std::allocator_traits<Alloc>::allocate(*this, 1);
  std::allocator_traits<Alloc>::construct(*this, new_entry,
                                          std::forward<Args>(args)...);

  for (;;) {
    for (bool pushing = false; not this->pushing.compare_exchange_weak(
             pushing, true, std::memory_order_acq_rel);)
      pushing = false;

    if (_size.load(std::memory_order_acquire) >= _max_size) {
      std::unique_lock lock(mxFull);

      if (_size.load(std::memory_order_acquire) >= _max_size) {
        pushing.store(false, std::memory_order_release);
        cvPop.wait(lock);
      } else
        pushing.store(false, std::memory_order_release);
    } else
      break;
  }

  for (bool popping = false; not this->popping.compare_exchange_weak(
           popping, true, std::memory_order_acq_rel);)
    popping = false;

  {
    entry_type* null_entry = nullptr;
    tail.compare_exchange_strong(null_entry, new_entry,
                                 std::memory_order_acq_rel);
  }

  auto current_head = head;
  new_entry->prev_entry = nullptr;
  new_entry->next_entry = current_head;

  if (current_head != nullptr)
    current_head->prev_entry = new_entry;

  head = new_entry;
  _size.fetch_add(1, std::memory_order_acq_rel);

  cvPush.notify_one();

  popping.store(false, std::memory_order_seq_cst);
  pushing.store(false, std::memory_order_seq_cst);
}

template <typename T, typename Alloc>
std::optional<T>
blocking_queue<T, Alloc>::pop() noexcept(
    std::is_nothrow_move_constructible_v<T>and
        std::is_nothrow_destructible_v<T>) {

  for (;;) {
    for (bool popping = false; not this->popping.compare_exchange_weak(
             popping, true, std::memory_order_acq_rel);)
      popping = false;

    entry_type* current_tail = tail.load(std::memory_order_acquire);
    if (current_tail == nullptr) {
      std::unique_lock lock(mxEmpty);

      if (tail.load(std::memory_order_acquire) == nullptr) {
        auto finished = _finished.load(std::memory_order_acquire);

        popping.store(false, std::memory_order_release);
        if (finished)
          return std::nullopt;
        else
          cvPush.wait(lock);
      } else
        popping.store(false, std::memory_order_release);
    } else
      break;
  }

  auto popped_entry = tail.exchange(nullptr, std::memory_order_acq_rel);

  assert(popped_entry->next_entry == nullptr);
  assert(_size.load(std::memory_order_acquire) > 0);
  _size.fetch_sub(1, std::memory_order_acq_rel);

  if (auto prev_entry = popped_entry->prev_entry; prev_entry == nullptr) {
    assert(head == popped_entry);

    head = nullptr;
  } else {
    prev_entry->next_entry = nullptr;
    tail.store(prev_entry, std::memory_order_release);
  }

  cvPop.notify_one();
  popping.store(false, std::memory_order_seq_cst);

  std::optional<T> out(std::in_place, std::move(popped_entry->entry));

#ifndef NDEBUG
  popped_entry->prev_entry = nullptr;
#endif
  std::allocator_traits<Alloc>::destroy(*this, popped_entry);
  std::allocator_traits<Alloc>::deallocate(*this, popped_entry, 1);

  return out;
}

template <typename T, typename Alloc>
std::optional<T>
blocking_queue<T, Alloc>::try_pop() noexcept(
    std::is_nothrow_move_constructible_v<T>and
        std::is_nothrow_destructible_v<T>) {

  for (bool popping = false; not this->popping.compare_exchange_weak(
           popping, true, std::memory_order_acq_rel);)
    popping = false;

  entry_type* current_tail = tail.load(std::memory_order_acquire);
  if (current_tail == nullptr) {
    popping.store(false, std::memory_order_release);
    return std::nullopt;
  }

  auto popped_entry = tail.exchange(nullptr, std::memory_order_acq_rel);

  assert(popped_entry->next_entry == nullptr);
  assert(_size.load(std::memory_order_acquire) > 0);
  _size.fetch_sub(1, std::memory_order_acq_rel);

  if (auto prev_entry = popped_entry->prev_entry; prev_entry == nullptr) {
    assert(head == popped_entry);

    head = nullptr;
  } else {
    prev_entry->next_entry = nullptr;
    tail.store(prev_entry, std::memory_order_release);
  }

  cvPop.notify_one();
  popping.store(false, std::memory_order_seq_cst);

  std::optional<T> out(std::in_place, std::move(popped_entry->entry));

#ifndef NDEBUG
  popped_entry->prev_entry = nullptr;
#endif
  std::allocator_traits<Alloc>::destroy(*this, popped_entry);
  std::allocator_traits<Alloc>::deallocate(*this, popped_entry, 1);

  return out;
}

} // namespace parallel

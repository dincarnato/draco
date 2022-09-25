#pragma once

#include "defaults.hpp"

#if USE_TBB
#include <tbb/parallel_for.h>

namespace parallel {

template <typename... Args> inline decltype(auto) parallel_for(Args &&...args) {
  return tbb::parallel_for(std::forward<Args>(args)...);
}

} /* namespace parallel */

#else

#include "blocked_range.hpp"

#include <thread>
#include <type_traits>
#include <vector>

namespace parallel {

template <typename Index, typename Func>
Func parallel_for(Index first, Index last, Index step, Func &&f) {
  std::vector<std::thread> threads;
  if (threadsPerLoop > 1) {
    assert(step >= 0);
    assert(step <= std::numeric_limits<unsigned>::max());
    auto subRanges =
        blocked_range<Index>(first, last)
            .split(threadsPerLoop - 1, static_cast<unsigned>(step));
    threads.reserve(subRanges.size());

    auto runner = [&]() -> decltype(auto) {
      if constexpr (std::is_invocable_v<Func, blocked_range<Index>>)
        return static_cast<Func &>(f);
      else {
        static_assert(std::is_invocable_v<Func, Index>);
        return [&f, step](blocked_range<Index> range) noexcept {
          for (auto index = std::begin(range); index < std::end(range);
               index += step)
            f(index);
        };
      }
    }();

    for (auto subRangeIter = std::begin(subRanges);
         subRangeIter < std::prev(std::end(subRanges)); ++subRangeIter)
      threads.emplace_back(runner, *subRangeIter);

    runner(subRanges.back());

    for (auto &thread : threads)
      thread.join();
  } else {
    if constexpr (std::is_invocable_v<Func, blocked_range<Index>>)
      f(blocked_range<Index>(first, last));
    else {
      static_assert(std::is_invocable_v<Func, Index>);
      for (Index index = first; index < last; index += step)
        f(index);
    }
  }

  return f;
}

template <typename Index, typename Func>
Func parallel_for(Index first, Index last, Func &&f) {
  return parallel_for(first, last, static_cast<Index>(1),
                      std::forward<Func>(f));
}

} /* namespace parallel */

#endif /* USE_TBB */

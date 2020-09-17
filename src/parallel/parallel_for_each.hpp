#pragma once

#include "defaults.hpp"

#if USE_TBB
#include <tbb/parallel_for_each.h>

namespace parallel {

template <typename... Args>
inline decltype(auto)
parallel_for_each(Args&&... args) {
  return tbb::parallel_for_each(std::forward<Args>(args)...);
}

} /* namespace parallel */

#else

#include "blocked_range.hpp"

#include <thread>
#include <vector>

namespace parallel {

template <typename Iterable, typename Func>
Func
parallel_for_each(Iterable&& iterable, Func&& f) {
  using diff_type = typename std::decay_t<Iterable>::difference_type;

  std::vector<std::thread> threads;

  if (threadsPerLoop > 1) {
    auto subRanges =
        blocked_range<diff_type>(
            0, std::distance(std::begin(iterable), std::end(iterable)))
            .split(threadsPerLoop - 1);
    threads.reserve(subRanges.size());

    auto runner = [&](blocked_range<diff_type> range) noexcept {
      auto&& iterableIter = std::next(std::begin(iterable), range[0]);
      auto&& end = std::end(range);
      for (auto index = std::begin(range); index < end; ++index, ++iterableIter)
        f(*iterableIter);
    };

    {
      auto&& subRangesEnd = std::prev(std::end(subRanges));
      for (auto subRangeIter = std::begin(subRanges);
           subRangeIter < subRangesEnd; ++subRangeIter)
        threads.emplace_back(runner, *subRangeIter);
    }

    runner(subRanges.back());

    for (auto& thread : threads)
      thread.join();
  } else {
    for (auto&& element : iterable)
      f(element);
  }

  return f;
}

} /* namespace parallel */

#endif /* USE_TBB */

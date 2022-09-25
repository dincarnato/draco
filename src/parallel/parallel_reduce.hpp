#pragma once

#include "defaults.hpp"

#if USE_TBB

#include <tbb/parallel_reduce.h>

namespace parallel {

template <typename... Args>
inline decltype(auto) parallel_reduce(Args &&...args) {
  return tbb::parallel_reduce(std::forward<Args>(args)...);
}

} /* namespace parallel */

#else

#include "blocked_range.hpp"

#include <future>
#include <vector>

namespace parallel {

template <typename Range, typename Value, typename Func, typename Reduction>
Value parallel_reduce(const Range &range, const Value &identity, Func func,
                      const Reduction &reduction) {
  std::vector<std::future<Value>> futures;
  auto subRanges = range.split(threadsPerLoop, 1);
  futures.reserve(subRanges.size());

  for (auto subRangeIter = std::begin(subRanges);
       subRangeIter < std::prev(std::end(subRanges)); ++subRangeIter)
    futures.emplace_back(
        std::async(std::launch::async, func, *subRangeIter, identity));

  Value returnValue = func(subRanges.back(), identity);
  for (auto &&future : futures)
    returnValue = reduction(future.get(), returnValue);

  return returnValue;
}

} /* namespace parallel */

#endif /* USE_TBB */

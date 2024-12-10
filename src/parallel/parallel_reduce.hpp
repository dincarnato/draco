#pragma once

#include "defaults.hpp"
#include <type_traits>

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
Value parallel_reduce(const Range &range, Value &&identity, Func func,
                      const Reduction &reduction) {
  using value_t = std::remove_cvref_t<Value>;

  std::vector<std::future<value_t>> futures;
  auto subRanges = range.split(threadsPerLoop, 1);
  futures.reserve(subRanges.size());

  auto &&local_identity = std::forward<Value>(identity);
  for (auto subRangeIter = std::begin(subRanges);
       subRangeIter < std::prev(std::end(subRanges)); ++subRangeIter) {
    if constexpr (std::is_invocable_v<Func, const Range &, value_t &&>) {
      futures.emplace_back(std::async(std::launch::async, func, *subRangeIter,
                                      value_t{local_identity}));
    } else {
      futures.emplace_back(
          std::async(std::launch::async, func, *subRangeIter, local_identity));
    }
  }

  value_t returnValue = ([&]() {
    if constexpr (std::is_invocable_v<Func, const Range &, value_t &&>) {
      return func(subRanges.back(), std::move(local_identity));
    } else {
      return func(subRanges.back(), local_identity);
    }
  })();
  for (auto &&future : futures) {
    if constexpr (std::is_invocable_v<Reduction, value_t &&, value_t &&>) {
      returnValue = reduction(future.get(), std::move(returnValue));
    } else {
      returnValue = reduction(future.get(), returnValue);
    }
  }

  return returnValue;
}

} /* namespace parallel */

#endif /* USE_TBB */

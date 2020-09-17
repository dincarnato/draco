#pragma once

#include "defaults.hpp"

#if USE_TBB

#include <tbb/blocked_range.h>

namespace parallel {

template <typename T>
using blocked_range = tbb::blocked_range<T>;

} /* namespace parallel */

#else

#include <vector>

namespace parallel {

template <typename T>
class blocked_range {
public:
  using value_type = T;
  using const_iterator = T;

  blocked_range(T first, T last) : first(first), last(last) {}

  const_iterator
  begin() const {
    return first;
  }

  const_iterator
  end() const {
    return last;
  }

  T
  size() const {
    return last - first;
  }

  std::vector<blocked_range>
  split(unsigned divisor, unsigned step = 1) const {
    std::vector<blocked_range> out;
    out.reserve(divisor);
    T divisorStep = size() / divisor;
    for (unsigned index = 0; index < divisor - 1; ++index) {
      T offset = index * divisorStep;
      T shift = offset % step;
      if (shift != 0)
        offset += step - shift;
      out.emplace_back(first + offset, first + (index + 1) * divisorStep);
    }

    T offset = (divisor - 1) * divisorStep;
    T shift = offset % step;
    if (shift != 0)
      offset += step - shift;
    out.emplace_back(first + offset, last);

    return out;
  }

  T operator[](std::size_t index) const { return first + index; }

private:
  T first;
  T last;
};

} /* namespace parallel */

#endif /* USE_TBB */

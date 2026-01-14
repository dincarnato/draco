#pragma once

#include "../args.hpp"

namespace test {

struct Args : ::Args {
  auto window_size() noexcept -> decltype(auto) { return (_window_size); }
  auto window_shift() noexcept -> decltype(auto) { return (_window_shift); }
  auto window_size_fraction_transcript_size() noexcept -> decltype(auto) {
    return (_window_size_fraction_transcript_size);
  }
};

} // namespace test

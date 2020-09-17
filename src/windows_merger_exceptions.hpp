#pragma once

#include <stdexcept>

namespace windows_merger::exceptions {

struct InvalidWindowSize : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct InvalidClustersSize : std::runtime_error {
  using std::runtime_error::runtime_error;
};

} // namespace windows_merger::exceptions

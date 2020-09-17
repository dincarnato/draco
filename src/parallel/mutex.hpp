#pragma once

#include "defaults.hpp"

#if USE_TBB

#include <tbb/mutex.h>

namespace parallel {
using mutex = tbb::mutex;
using lock_guard = tbb::mutex::scoped_lock;
} /* namespace parallel */

#else

#include <mutex>
namespace parallel {
using mutex = std::mutex;
using lock_guard = std::lock_guard<std::mutex>;
} /* namespace parallel */

#endif /* USE_TBB */

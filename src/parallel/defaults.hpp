#pragma once

#ifndef USE_TBB
#define USE_TBB 0
#endif

#ifndef SINGLE_THREAD
#define SINGLE_THREAD 0
#endif

#if not USE_TBB

namespace parallel {
#if SINGLE_THREAD
constexpr unsigned threadsPerLoop = 1;
#else
constexpr unsigned threadsPerLoop = 4;
#endif
} /* namespace parallel */

#endif

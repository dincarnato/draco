#pragma once

#include <atomic>
#include <concepts>
#include <exception>
#include <fmt/base.h>
#include <fmt/color.h>
#include <iterator>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>

#include "fmt/base.h"
#include "fmt/color.h"

namespace logger {

enum class Level {
  trace,
  debug,
  info,
  warn,
  error,
};

constexpr std::optional<Level>
parse_level(std::string_view log_level_str) noexcept {
  if (log_level_str == "trace") {
    return std::optional{logger::Level::trace};
  } else if (log_level_str == "debug") {
    return std::optional{logger::Level::debug};
  } else if (log_level_str == "info") {
    return std::optional{logger::Level::info};
  } else if (log_level_str == "warn") {
    return std::optional{logger::Level::warn};
  } else if (log_level_str == "error") {
    return std::optional{logger::Level::error};
  } else {
    return std::optional<logger::Level>{};
  }
}

struct Logger {
  constexpr Logger() = default;
  constexpr Logger(Level level) noexcept : level(level) {}

  void set_level(Level level) noexcept {
    this->level.store(level, std::memory_order::relaxed);
  }

  Level get_level() const noexcept {
    return this->level.load(std::memory_order::relaxed);
  }

private:
  std::atomic<Level> level{Level::info};
};

extern Logger instance;

template <typename... Args>
inline void log(Level level, fmt::format_string<Args...> format_string,
                Args &&...args) noexcept {
  if (level < instance.get_level()) {
    return;
  }

  auto [text_style, level_str] = ([=]() {
    switch (level) {
    case Level::trace:
      return std::pair{fmt::fg(fmt::color::gray), "TRACE"};
    case Level::debug:
      return std::pair{fmt::fg(fmt::color::white), "DEBUG"};
    case Level::info:
      return std::pair{fmt::fg(fmt::color::cyan), "INFO"};
    case Level::warn:
      return std::pair{fmt::emphasis::bold | fmt::fg(fmt::color::orange),
                       "WARN"};
    case Level::error:
      return std::pair{fmt::emphasis::bold | fmt::fg(fmt::color::red), "ERROR"};
    }

    // Unreachable
    std::terminate();
  })();

  std::string s =
      fmt::format("[{}]: ", fmt::styled(level_str, std::move(text_style)));
  fmt::format_to(std::back_inserter(s), format_string,
                 std::forward<Args>(args)...);
  fmt::println(stderr, "{}", s);
}

template <typename... Args, std::invocable<Args...> Fn>
inline auto on_log_level(Level level, Fn &&fn, Args &&...args)
    -> decltype(auto) {
  using return_t = std::invoke_result_t<Fn, Args...>;

  auto run_fn = level >= instance.get_level();

  if constexpr (std::is_void_v<return_t>) {
    if (run_fn) {
      fn(std::forward<Args>(args)...);
    }
  } else {
    if (run_fn) {
      return std::optional(fn(std::forward<Args>(args)...));
    } else {
      return std::optional<return_t>();
    }
  }
}

#define IMPL_LOG_FUNS_WITH_LEVEL(LEVEL)                                        \
  template <typename... Args>                                                  \
  inline void LEVEL(fmt::format_string<Args...> format_string,                 \
                    Args &&...args) noexcept {                                 \
    log(Level::LEVEL, std::move(format_string), std::forward<Args>(args)...);  \
  }                                                                            \
  template <typename... Args, std::invocable<Args...> Fn>                      \
  inline auto on_##LEVEL##_level(Fn &&fn, Args &&...args) -> decltype(auto) {  \
    on_log_level(Level::LEVEL, std::forward<Fn>(fn),                           \
                 std::forward<Args>(args)...);                                 \
  }

IMPL_LOG_FUNS_WITH_LEVEL(trace);
IMPL_LOG_FUNS_WITH_LEVEL(debug);
IMPL_LOG_FUNS_WITH_LEVEL(info);
IMPL_LOG_FUNS_WITH_LEVEL(warn);
IMPL_LOG_FUNS_WITH_LEVEL(error);

} // namespace logger

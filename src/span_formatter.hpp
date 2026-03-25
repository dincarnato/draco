#pragma once

#include "fmt/base.h"

#include <ranges>
#include <span>
#include <vector>

template <typename T> struct SpanFormatter {
  constexpr SpanFormatter() = default;
  constexpr explicit SpanFormatter(std::span<const T> span) noexcept
      : data_(span) {}
  constexpr explicit SpanFormatter(std::span<T> span) noexcept : data_(span) {}
  constexpr explicit SpanFormatter(std::vector<T> span) noexcept
      : data_(span) {}

  constexpr std::span<const T> data() const noexcept { return data_; }

protected:
  std::span<const T> data_;
};

template <typename T> struct fmt::formatter<SpanFormatter<T>> {
  fmt::format_context::iterator format(const SpanFormatter<T> &f,
                                       fmt::format_context &ctx) const {
    auto out_iter = ctx.out();

    ctx.advance_to(fmt::format_to(out_iter, "["));
    auto data = f.data();
    auto iter = std::ranges::begin(data);
    if (iter != std::ranges::end(data)) {
      ctx.advance_to(fmt::format_to(out_iter, "{}", *iter));

      for (auto const &element :
           std::ranges::subrange(std::next(iter), std::ranges::end(data))) {
        ctx.advance_to(fmt::format_to(out_iter, ", {}", element));
      }
    }

    return fmt::format_to(out_iter, "]");
  }

  constexpr auto parse(fmt::format_parse_context &ctx) { return ctx.begin(); }
};

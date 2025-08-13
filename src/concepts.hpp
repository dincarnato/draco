#pragma once

#include <concepts>
#include <utility>

template <typename Fn, typename R, typename... Args>
concept InvocableR = requires(Fn &&fn, Args &&...args) {
  requires std::invocable<Fn, Args...>;
  { fn(std::forward<Args>(args)...) } -> std::same_as<R>;
};

static_assert(InvocableR<decltype([](unsigned, unsigned) { return 3u; }),
                         unsigned, unsigned, unsigned>);
static_assert(!InvocableR<decltype([](unsigned, unsigned) { return 3u; }), int,
                          unsigned, unsigned>);
static_assert(!InvocableR<decltype([](unsigned, unsigned) { return 3u; }), void,
                          unsigned, unsigned>);
static_assert(!InvocableR<decltype([](unsigned) { return 3u; }), unsigned,
                          unsigned, unsigned>);

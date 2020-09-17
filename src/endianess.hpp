#pragma once

#include <cstdint>
#include <type_traits>

enum class Endianess { little, big };

struct system_is_little_endian {
  static constexpr bool value = (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__);
};

struct system_is_big_endian {
  static constexpr bool value = (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__);
};

namespace detail {

template <std::size_t index, class T>
inline constexpr std::enable_if_t<index == 0, T>
swapBytes(const T&) {
  return 0;
}

template <std::size_t index, class T>
inline constexpr std::enable_if_t<index != 0, T>
swapBytes(const T& t) {
  return swapBytes<index - 1>(t) |
         (((t >> ((sizeof(T) - index) * 8)) & 0xff) << ((index - 1) * 8)) |
         (((t >> ((index - 1) * 8)) & 0xff) << ((sizeof(T) - index) * 8));
}

} /* namespace detail */

template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
constexpr void
swapBytes(T& t) {
  t = detail::swapBytes<sizeof(T) / 2>(t);
}

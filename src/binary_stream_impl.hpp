#pragma once

#include "binary_stream.hpp"

#include <iomanip>

template <typename Stream>
BinaryStream<Stream>::BinaryStream(Stream& stream) : stream(&stream) {}

template <typename Stream>
template <typename T, typename Out>
std::enable_if_t<system_is_little_endian::value, Out>
BinaryStream<Stream>::operator>>(T& t) {
  stream->read(reinterpret_cast<char*>(&t), sizeof(T));
  return *this;
}

template <typename Stream>
template <typename T, typename Out>
std::enable_if_t<system_is_big_endian::value, Out>
BinaryStream<Stream>::operator>>(T& t) {
  stream->read(reinterpret_cast<char*>(&t), sizeof(T));
  swapBytes(t);
  return *this;
}

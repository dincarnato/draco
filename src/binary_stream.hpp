#pragma once

#include "endianess.hpp"

template <typename Stream>
class BinaryStream {
public:
  BinaryStream(Stream& stream);

  template <typename T, typename Out = BinaryStream&>
  std::enable_if_t<system_is_little_endian::value, Out> operator>>(T& t);

  template <typename T, typename Out = BinaryStream&>
  std::enable_if_t<system_is_big_endian::value, Out> operator>>(T& t);

private:
  Stream* stream;
};

#include "binary_stream_impl.hpp"

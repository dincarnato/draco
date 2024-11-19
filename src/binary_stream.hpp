#pragma once

#include "endianess.hpp"
#include <ios>

template <typename Stream> class BinaryStream {
public:
  using pos_type = typename Stream::pos_type;
  using off_type = typename Stream::off_type;

  BinaryStream(Stream &stream);

  template <typename T, typename Out = BinaryStream &>
  std::enable_if_t<system_is_little_endian::value, Out> operator>>(T &t);

  template <typename T, typename Out = BinaryStream &>
  std::enable_if_t<system_is_big_endian::value, Out> operator>>(T &t);

  template <typename T, typename Out = BinaryStream &>
  std::enable_if_t<system_is_little_endian::value, Out> operator<<(T const &t);

  template <typename T, typename Out = BinaryStream &>
  std::enable_if_t<system_is_big_endian::value, Out> operator<<(T t);

  bool is_open() const;
  pos_type tellp();
  BinaryStream<Stream> &seekp(pos_type pos);
  BinaryStream<Stream> &seekp(off_type off, std::ios_base::seekdir dir);

  Stream &inner();
  Stream const &inner() const;

private:
  Stream *stream;
};

#include "binary_stream_impl.hpp"

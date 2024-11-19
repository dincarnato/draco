#pragma once

#include "binary_stream.hpp"

#include <iomanip>

template <typename Stream>
BinaryStream<Stream>::BinaryStream(Stream &stream) : stream(&stream) {}

template <typename Stream>
template <typename T, typename Out>
std::enable_if_t<system_is_little_endian::value, Out>
BinaryStream<Stream>::operator>>(T &t) {
  stream->read(reinterpret_cast<char *>(&t), sizeof(T));
  return *this;
}

template <typename Stream>
template <typename T, typename Out>
std::enable_if_t<system_is_big_endian::value, Out>
BinaryStream<Stream>::operator>>(T &t) {
  stream->read(reinterpret_cast<char *>(&t), sizeof(T));
  swapBytes(t);
  return *this;
}

template <typename Stream>
template <typename T, typename Out>
std::enable_if_t<system_is_little_endian::value, Out>
BinaryStream<Stream>::operator<<(T const &t) {
  stream->write(reinterpret_cast<char const *>(&t), sizeof(T));
  return *this;
}

template <typename Stream>
template <typename T, typename Out>
std::enable_if_t<system_is_big_endian::value, Out>
BinaryStream<Stream>::operator<<(T t) {
  swapBytes(t);
  stream->write(reinterpret_cast<char const *>(&t), sizeof(T));
  return *this;
}

template <typename Stream> bool BinaryStream<Stream>::is_open() const {
  return stream->is_open();
}

template <typename Stream> auto BinaryStream<Stream>::tellp() -> pos_type {
  return stream->tellp();
}

template <typename Stream>
BinaryStream<Stream> &BinaryStream<Stream>::seekp(pos_type pos) {
  stream->seekp(pos);
  return *this;
}

template <typename Stream>
BinaryStream<Stream> &BinaryStream<Stream>::seekp(off_type off,
                                                  std::ios_base::seekdir dir) {
  stream->seekp(off, dir);
  return *this;
}

template <typename Stream> Stream &BinaryStream<Stream>::inner() {
  return *stream;
}

template <typename Stream> Stream const &BinaryStream<Stream>::inner() const {
  return *stream;
}

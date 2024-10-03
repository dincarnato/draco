#pragma once

#include "tokenizer_iterator.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <type_traits>
#include <utility>

#if __cpp_lib_to_chars
#include <charconv>
#endif

template <typename String>
TokenizerIterator<String>::TokenizerIterator(const String &str, char token)
    : str(&str), token(token), iter(&*std::begin(str)), end(&*std::end(str)),
      next((iter < end) ? std::find(iter, end, token) : end) {}

template <typename String> TokenizerIterator<String>::operator bool() const {
  return iter < end;
}

template <typename String>
typename TokenizerIterator<String>::out_type
TokenizerIterator<String>::operator*() const {
  return out_type(iter, static_cast<typename out_type::size_type>(next - iter));
}

template <typename String>
typename TokenizerIterator<String>::out_type *
TokenizerIterator<String>::operator->() const {
  static out_type out;
  out = out_type(iter, static_cast<typename out_type::size_type>(next - iter));
  return &out;
}

template <typename String>
TokenizerIterator<String> &TokenizerIterator<String>::operator++() {
  iter = std::next(next);
  if (iter != end)
    next = std::find(iter, end, token);

  return *this;
}

template <typename String>
TokenizerIterator<String> TokenizerIterator<String>::operator++(int) {
  TokenizerIterator copy(*this);
  this->operator++();
  return copy;
}

template <typename String> int TokenizerIterator<String>::toi() const {
#if __cpp_lib_to_chars
  int out;
  auto result = std::from_chars(iter, next, out);
  if (result.ptr == iter || result.ec == std::errc::result_out_of_range) {
    throw std::runtime_error("unable to convert from characters to int");
  }
  return out;
#else
  return std::stoi(std::string(iter, next));
#endif
}

template <typename String> unsigned TokenizerIterator<String>::tou() const {
#if __cpp_lib_to_chars
  unsigned out;
  auto result = std::from_chars(iter, next, out);
  if (result.ptr == iter || result.ec == std::errc::result_out_of_range) {
    throw std::runtime_error(
        "unable to convert from characters to unsigned int");
  }
  return out;
#else
  return static_cast<unsigned>(std::stoi(std::string(iter, next)));
#endif
}

template <typename String>
TokenizerIterator<std::decay_t<String>> makeTokenizerIterator(String &&str,
                                                              char token) {
  return TokenizerIterator<std::decay_t<String>>(std::forward<String>(str),
                                                 token);
}

template <typename String>
TokenizerIterator<String> &
TokenizerIterator<String>::operator+=(std::ptrdiff_t offset) {
  assert(offset >= 0);
  for (std::ptrdiff_t index = 0; index < offset; ++index)
    operator++();
  return *this;
}

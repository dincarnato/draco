#pragma once

#if __cplusplus > 201402L
#if __has_include(<string_view>)
#include <string_view>
namespace {
using string_out_type = std::string_view;
}
#elif __has_include(<experimental/string_view>)
#include <experimental/string_view>
namespace {
using string_out_type = std::experimental::string_view;
}
#else
#include <string>
namespace {
using string_out_type = std::string;
}
#endif
#else
#include <string>
namespace {
using string_out_type = std::string;
}
#endif

template <typename String> class TokenizerIterator {
public:
  using string_type = String;
  using string_pointer_type = typename String::const_pointer;
  using out_type = string_out_type;

  explicit TokenizerIterator(const String &str, char token = ' ');
  TokenizerIterator &operator++();
  TokenizerIterator operator++(int);
  TokenizerIterator &operator+=(std::ptrdiff_t offset);
  out_type operator*() const;
  out_type *operator->() const;
  int toi() const;
  unsigned tou() const;

  operator bool() const;

private:
  const String *str;
  char token;
  string_pointer_type iter, end, next;
};

template <typename String>
TokenizerIterator<std::decay_t<String>> makeTokenizerIterator(String &&str,
                                                              char token = ' ');

#include "tokenizer_iterator_impl.hpp"

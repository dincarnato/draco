#pragma once

#include "args.hpp"
#include "math.hpp"

#include <array>
#include <functional>
#include <type_traits>
#include <vector>

#include <armadillo>

class RingmapData;

namespace ringmap {
using umat = ::arma::Mat<unsigned>;
using uvec = ::arma::Col<unsigned>;
} /* namespace ringmap */

std::vector<std::vector<unsigned>>
cleanStructures(const RingmapData &ringmapData,
                const std::vector<std::vector<unsigned>> &assignments,
                Args const &args);

constexpr std::size_t operator""_sz(unsigned long long n) { return n; }

namespace detail {

template <unsigned k, unsigned depth, typename Iterable>
inline std::enable_if_t<
    depth == k - 1,
    typename std::vector<
        std::array<typename std::decay_t<Iterable>::value_type, k>>::iterator>
combinations(
    Iterable &&elements,
    typename std::vector<std::array<typename std::decay_t<Iterable>::value_type,
                                    k>>::iterator &&outIter,
    unsigned elementIndex) {
  for (auto elementIter =
                std::ranges::next(std::ranges::begin(elements), elementIndex),
            endIter = std::ranges::end(elements);
       elementIter < endIter; ++elementIter, ++outIter)
    (*outIter)[depth] = *elementIter;
  return outIter;
}

template <unsigned k, unsigned depth, typename Iterable>
inline std::enable_if_t<
    (depth < k - 1),
    typename std::vector<
        std::array<typename std::decay_t<Iterable>::value_type, k>>::iterator>
combinations(
    Iterable &&elements,
    typename std::vector<std::array<typename std::decay_t<Iterable>::value_type,
                                    k>>::iterator &&outIter,
    unsigned elementIndex) {
  for (unsigned index = elementIndex; index < elements.size(); ++index) {
    const auto &element = elements[index];
    std::size_t elementRepetitions =
        binomial(elements.size() - index - 1, k - depth - 1);

    auto &&currentEnd = std::ranges::next(outIter, elementRepetitions);
    for (auto iter = outIter; iter < currentEnd; ++iter)
      (*iter)[depth] = element;

    {
      auto iter = outIter;
      do {
        iter = combinations<k, depth + 1>(elements, std::move(iter), index + 1);
      } while (iter < currentEnd);
    }

    outIter = std::move(currentEnd);
  }

  return outIter;
}

} /* namespace detail */

template <unsigned k, typename Iterable>
auto combinations(Iterable &&elements)
    -> std::vector<std::array<typename std::decay_t<Iterable>::value_type, k>> {
  using value_type = typename std::decay_t<Iterable>::value_type;
  std::vector<std::array<value_type, k>> out(binomial(elements.size(), k));
  detail::combinations<k, 0u>(std::forward<Iterable>(elements),
                              std::ranges::begin(out), 0);
  return out;
}

template <typename InputIt1, typename InputIt2, typename Compare>
std::size_t count_intersections(InputIt1 first1, InputIt1 last1,
                                InputIt2 first2, InputIt2 last2, Compare comp) {
  if (first1 == last1 or first2 == last2)
    return 0;

  std::size_t count = 0;
  for (;;) {
    first1 = std::ranges::lower_bound(first1, last1, *first2, comp);
    if (first1 == last1)
      break;

    first2 = std::ranges::lower_bound(first2, last2, *first1, comp);
    if (first2 == last2)
      break;

    if (*first1 == *first2) {
      ++count;
      if (++first1 == last1 or ++first2 == last2)
        break;
    }
  }

  return count;
}

namespace arma {

template <typename T> inline bool operator==(const Row<T> &a, const Row<T> &b) {
  auto size = a.size();
  assert(size == b.size());
  for (decltype(size) index = 0; index < size; ++index) {
    if (a[index] != b[index])
      return false;
  }
  return true;
}

template <typename T> inline bool operator==(const Col<T> &a, const Col<T> &b) {
  auto size = a.size();
  assert(size == b.size());
  for (decltype(size) index = 0; index < size; ++index) {
    if (a[index] != b[index])
      return false;
  }
  return true;
}

template <typename T> inline bool operator!=(const Row<T> &a, const Row<T> &b) {
  return not(a == b);
}

template <typename T> inline bool operator!=(const Col<T> &a, const Col<T> &b) {
  return not(a == b);
}

} /* namespace arma */

template <typename InputIt1, typename InputIt2>
inline std::size_t count_intersections(InputIt1 first1, InputIt1 last1,
                                       InputIt2 first2, InputIt2 last2) {
  return count_intersections(std::move(first1), std::move(last1),
                             std::move(first2), std::move(last2),
                             std::less<typename InputIt1::value_type>());
}

/*
 * Helper functions for debugging armadillo matrices and vectors
 */
#ifndef NDEBUG

template <typename T>
T arma_get_vec(const arma::Col<T> &vec, std::size_t index) {
  return vec[index];
}

template <typename T>
T arma_get_rowvec(const arma::Row<T> &vec, std::size_t index) {
  return vec[index];
}

template <typename T>
T arma_get_mat(const arma::Mat<T> &mat, std::size_t row, std::size_t col) {
  return mat(row, col);
}

#define DEFINE_ARMA_UTIL(type)                                                 \
  template type arma_get_vec<type>(const arma::Col<type> &, std::size_t);      \
  template type arma_get_rowvec<type>(const arma::Row<type> &, std::size_t);   \
  template type arma_get_mat<type>(const arma::Mat<type> &, std::size_t,       \
                                   std::size_t);

DEFINE_ARMA_UTIL(double)
DEFINE_ARMA_UTIL(float)
DEFINE_ARMA_UTIL(long)
DEFINE_ARMA_UTIL(unsigned long)
DEFINE_ARMA_UTIL(int)
DEFINE_ARMA_UTIL(unsigned)

#endif

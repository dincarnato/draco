#pragma once

#include "nostd/type_traits.hpp"

#include <armadillo>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <memory>
#include <numeric>
#include <range/v3/iterator/concepts.hpp>
#include <type_traits>
#include <vector>

template <typename Iterable> constexpr double mean(Iterable &&iterable);

template <typename Iterable> constexpr double variance(Iterable &&iterable);

template <typename Iterable> constexpr double stddev(Iterable &&iterable);

template <typename Iterable> std::vector<float> ranks(Iterable &&iterable);

template <typename IterableA, typename IterableB>
double covariance(IterableA &&iterableA, IterableB &&iterableB);

enum class CorrelationType { spearman };

template <CorrelationType correlationType = CorrelationType::spearman,
          typename IterableA, typename IterableB>
auto correlation(IterableA &&iterableA, IterableB &&iterableB)
    -> std::enable_if_t<correlationType == CorrelationType::spearman, double>;

constexpr std::uint64_t factorial(int n) noexcept(false) {
  if (n < 0)
    throw std::logic_error("n must be positive");

  if (n > 22)
    throw std::overflow_error(
        "max factoriable number without overflowing 64bit: 22");

  std::int64_t out = 1;
  for (int index = 2; index <= n; ++index)
    out *= index;
  return static_cast<std::uint64_t>(out);
}

namespace detail {
template <typename T, typename = void> struct numeric_traits;

template <typename T>
struct numeric_traits<T, std::enable_if_t<std::is_integral_v<T>>> {
  using signed_type = std::make_signed_t<T>;
  using unsigned_type = std::make_unsigned_t<T>;
  using signed_integral_type = std::make_signed_t<T>;
};

template <typename T>
struct numeric_traits<T, std::enable_if_t<std::is_floating_point_v<T>>> {
  using signed_type = T;
  using unsigned_type = T;
  using signed_integral_type = std::int64_t;
};

} // namespace detail

template <typename T, std::ptrdiff_t max_num_factors = 64>
constexpr typename detail::numeric_traits<T>::unsigned_type
binomial_generic(unsigned n, unsigned k) {
  if (k == 0 or k >= n)
    return 1;

  using result_t = typename detail::numeric_traits<T>::signed_type;
  using working_t = typename detail::numeric_traits<T>::signed_integral_type;

  const std::ptrdiff_t n_factors = n - k;
  if (n_factors > max_num_factors)
    throw std::logic_error("n - k > max_num_factors; check your code or "
                           "consider increasing max_num_factors.");

  std::array<working_t, static_cast<std::size_t>(max_num_factors)>
      num_factors{};
  const auto end_num_factors = std::next(std::begin(num_factors), n_factors);
  {
    working_t factor = k + 1;
    for (auto num_factor_iter = std::begin(num_factors);
         num_factor_iter != end_num_factors; ++num_factor_iter, ++factor)
      *num_factor_iter = factor;
  }

  for (working_t den_factor = 2; den_factor < (n - k + 1); ++den_factor) {
    auto factor = den_factor;
    for (auto num_factor_iter = std::begin(num_factors);
         num_factor_iter != end_num_factors; ++num_factor_iter) {
      auto &num_factor = *num_factor_iter;
      auto divider = std::gcd(num_factor, factor);
      if (divider != 0) {
        num_factor /= divider;
        factor /= divider;

        if (factor == 1)
          break;
      }
    }

    assert(factor == 1);
  }

  result_t binomial = 1;
  for (auto num_factor_iter = std::begin(num_factors);
       num_factor_iter != end_num_factors; ++num_factor_iter) {
    if constexpr (std::is_integral_v<result_t>) {
      auto new_binomial = binomial * static_cast<result_t>(*num_factor_iter);
      if (new_binomial < binomial)
        throw std::overflow_error("binomial overflows");
      binomial = new_binomial;
    } else {
      binomial *= static_cast<result_t>(*num_factor_iter);
      if (std::isinf(binomial))
        throw std::overflow_error("binomial reached infinity");
    }
  }
  return static_cast<typename detail::numeric_traits<T>::unsigned_type>(
      binomial);
}

constexpr std::uint64_t binomial(unsigned n, unsigned k) {
  return binomial_generic<std::uint64_t>(n, k);
}

constexpr double binomial_floating(unsigned n, unsigned k) {
  return binomial_generic<double, 1024>(n, k);
}

template <typename Iterable> constexpr double mean(Iterable &&iterable) {
  std::size_t size = iterable.size();
  return static_cast<double>(
             std::accumulate(std::begin(std::forward<Iterable>(iterable)),
                             std::end(std::forward<Iterable>(iterable)),
                             typename std::decay_t<Iterable>::value_type(0))) /
         size;
}

template <typename Iterable> constexpr double variance(Iterable &&iterable) {
  double meanValue = mean(iterable);
  std::size_t size = iterable.size();
  return static_cast<double>(std::accumulate(
             std::begin(std::forward<Iterable>(iterable)),
             std::end(std::forward<Iterable>(iterable)), 0.,
             [meanValue](double accumulator, auto &&value) {
               return accumulator + std::pow(value - meanValue, 2);
             })) /
         (size - 1);
}

template <typename Iterable> constexpr double stddev(Iterable &&iterable) {
  return std::sqrt(variance(std::forward<Iterable>(iterable)));
}

template <typename Iterable> std::vector<float> ranks(Iterable &&iterable) {
  // TODO: if constexpr to work with objects without operator[] and/or using
  // iterators/pointers for order
  std::vector<std::size_t> order(iterable.size());
  std::iota(std::begin(order), std::end(order), std::size_t(0));
  std::sort(std::begin(order), std::end(order),
            [&iterable](std::size_t a, std::size_t b) {
              return iterable[a] < iterable[b];
            });

  std::vector<float> ranks(order.size());
  std::size_t index = 0;
  for (auto indexIter = std::begin(order); indexIter < std::end(order);) {
    const auto &currentValue = iterable[*indexIter];
    auto otherIndexIter = std::next(indexIter);
    for (; otherIndexIter < std::end(order); ++otherIndexIter) {
      if (iterable[*otherIndexIter] != currentValue)
        break;
    }

    auto indices =
        static_cast<unsigned>(std::distance(indexIter, otherIndexIter));
    float rankIndex = static_cast<float>(index * 2 + indices - 1) / 2.f;
    for (; indexIter < otherIndexIter; ++indexIter)
      ranks[*indexIter] = rankIndex;
    index += indices;
  }

  return ranks;
}

template <typename IterableA, typename IterableB>
double covariance(IterableA &&iterableA, IterableB &&iterableB) {
  std::size_t size = iterableA.size();
  assert(size == iterableB.size());

  double meanA = mean(iterableA);
  double meanB = mean(iterableB);

  double covariance = 0.;

  {
    auto iterB = std::begin(iterableB);
    auto endA = std::end(iterableA);
    for (auto iterA = std::begin(iterableA); iterA != endA; ++iterA, ++iterB)
      covariance += (*iterA - meanA) * (*iterB - meanB);
  }
  return covariance / static_cast<double>(size - 1);
}

template <CorrelationType correlationType, typename IterableA,
          typename IterableB>
auto correlation(IterableA &&iterableA, IterableB &&iterableB)
    -> std::enable_if_t<correlationType == CorrelationType::spearman, double> {
  auto &&rankA = ranks(std::forward<IterableA>(iterableA));
  auto &&rankB = ranks(std::forward<IterableB>(iterableB));
  assert(rankA.size() > 1);
  assert(rankB.size() > 1);
  assert(std::any_of(
      std::next(std::begin(rankA)), std::end(rankA),
      [first = rankA.front()](double value) { return value != first; }));
  assert(std::any_of(
      std::next(std::begin(rankB)), std::end(rankB),
      [first = rankB.front()](double value) { return value != first; }));

  const double divisor = stddev(rankA) * stddev(rankB);
  assert(divisor != 0.);

  return covariance(std::move(rankA), std::move(rankB)) / divisor;
}

template <typename PredMat, typename TrueMat>
double matthewCorrelationCoefficient(PredMat &&predictionMatrix,
                                     TrueMat &&trueMatrix) {
  unsigned nRows = predictionMatrix.n_rows;

  assert(nRows == trueMatrix.n_rows);
  assert(predictionMatrix.n_cols == trueMatrix.n_cols);

  std::decay_t<PredMat> predictionDispersion =
      predictionMatrix - arma::repmat(arma::mean(predictionMatrix), nRows, 1);
  std::decay_t<TrueMat> trueDispersion =
      trueMatrix - arma::repmat(arma::mean(trueMatrix), nRows, 1);

  auto cov = [](const auto &matA, const auto &matB) {
    return arma::sum(arma::sum(matA % matB)) / matA.n_cols;
  };

  double covPrediction = cov(predictionDispersion, predictionDispersion);
  if (covPrediction == 0.)
    return -1.;

  double covTrue = cov(trueDispersion, trueDispersion);
  if (covTrue == 0.)
    return -1.;

  return cov(predictionDispersion, trueDispersion) /
         std::sqrt(covPrediction * covTrue);
}

template <typename Iterable>
auto getPercentile(Iterable &&iterable,
                   float percentile) -> std::decay_t<decltype(iterable[0])> {
  std::remove_reference_t<Iterable> *usableIterable = &iterable;
  std::unique_ptr<std::decay_t<Iterable>> localIterable;
  if (not std::is_sorted(std::begin(iterable), std::end(iterable))) {
    localIterable = std::make_unique<std::decay_t<Iterable>>(
        std::forward<Iterable>(iterable));
    std::sort(std::begin(*localIterable), std::end(*localIterable));
    usableIterable = localIterable.get();
  }

  std::size_t index = std::min(usableIterable->size() - 1,
                               static_cast<decltype(usableIterable->size())>(
                                   percentile / 100. * usableIterable->size()));
  assert(index < usableIterable->size());
  return (*usableIterable)[index];
}

template <typename Iterable>
std::decay_t<Iterable> getWinsorized(Iterable &&iterable, float percentile) {
  std::decay_t<Iterable> out(std::forward<Iterable>(iterable));
  if (out.size() == 0)
    return out;

  auto maxValue = getPercentile(out, percentile);
  std::transform(
      std::begin(out), std::end(out), std::begin(out),
      [maxValue](auto value) -> std::decay_t<decltype(*std::begin(out))> {
        if (value >= maxValue)
          return 1;
        return static_cast<double>(value) / maxValue;
      });

  return out;
}

template <typename IterA, typename IterB, typename ScoreFun,
          typename CompareFun>
std::vector<std::size_t>
getBestMatchingIndices(IterA refBegin, IterA refEnd, IterB matchBegin,
                       IterB matchEnd, ScoreFun &&scoreFun,
                       CompareFun &&compareFun) {
  using ref_reference_type = typename std::iterator_traits<IterA>::reference;
  using match_reference_type = typename std::iterator_traits<IterB>::reference;
#if defined(__cpp_lib_is_invocable) && __cpp_lib_is_invocable >= 201703
  using score_type =
      std::invoke_result_t<std::decay_t<ScoreFun>, ref_reference_type,
                           match_reference_type>;
#else
  using score_type = std::result_of_t<std::decay_t<ScoreFun>(
      ref_reference_type, match_reference_type)>;
#endif

  const auto refElements = static_cast<int>(std::distance(refBegin, refEnd));
  const auto matchingElements =
      static_cast<int>(std::distance(matchBegin, matchEnd));
  arma::Mat<score_type> scores(factorial(matchingElements),
                               static_cast<arma::uword>(refElements));
#ifndef NDEBUG
  std::fill(std::begin(scores), std::end(scores),
            std::numeric_limits<score_type>::max());
#endif

  std::vector<int> matchIndices(static_cast<std::size_t>(matchingElements));
  const auto matchIndicesEnd = std::end(matchIndices);
  std::iota(std::begin(matchIndices), std::end(matchIndices), 0);
  std::vector<std::vector<std::size_t>> allIndices(
      scores.n_rows, std::vector<std::size_t>(static_cast<std::size_t>(
                         std::min(refElements, matchingElements))));

  std::size_t matchesSize = 0;
  {
    unsigned long permutationIndex = 0;
    do {
      assert(permutationIndex <
             std::numeric_limits<decltype(permutationIndex)>::max());
      auto refIter = refBegin;
      auto matchIndexIter = std::begin(matchIndices);

      std::size_t elementIndex = 0;
      for (; refIter != refEnd and matchIndexIter != matchIndicesEnd;
           ++refIter, ++elementIndex, ++matchIndexIter) {
        assert(elementIndex < scores.n_cols);
        assert(permutationIndex < scores.n_rows);
        assert(*matchIndexIter < matchingElements);
        if constexpr (ranges::RandomAccessIterator<
                          std::remove_reference_t<IterA>>)
          assert(refIter < refEnd);
        scores(permutationIndex, elementIndex) = scoreFun(
            *refIter, *std::next(matchBegin,
                                 static_cast<std::ptrdiff_t>(*matchIndexIter)));
      }
      std::copy(std::begin(matchIndices), matchIndexIter,
                std::begin(allIndices[permutationIndex]));
      matchesSize = elementIndex;

      ++permutationIndex;
    } while (
        next_permutation(std::begin(matchIndices), std::end(matchIndices)));

    assert(permutationIndex == factorial(matchingElements));
  }

#ifndef NDEBUG
  for (std::size_t row = 0; row < scores.n_rows; ++row) {
    score_type accumulator = 0;
    for (std::size_t col = 0; col < matchesSize; ++col) {
      score_type score = scores(row, col);
      assert(score != std::numeric_limits<score_type>::max());

      if (score >= 0)
        assert(accumulator < std::numeric_limits<score_type>::max() - score);
      else
        assert(accumulator > std::numeric_limits<score_type>::min() - score);
      accumulator += score;
    }

    (void)accumulator;
  }
#endif

  const arma::Col<score_type> permutationsScore =
      arma::sum(scores.head_cols(matchesSize), 1);
  const auto bestIndex = static_cast<std::size_t>(
      std::distance(std::begin(permutationsScore),
                    std::max_element(std::begin(permutationsScore),
                                     std::end(permutationsScore),
                                     std::forward<CompareFun>(compareFun))));
  return allIndices[bestIndex];
}

template <typename IterA, typename IterB, typename ScoreFun>
std::vector<std::size_t>
getBestMatchingIndices(IterA refBegin, IterA refEnd, IterB matchBegin,
                       IterB matchEnd, ScoreFun &&scoreFun) {
  using ref_reference_type = typename std::iterator_traits<IterA>::reference;
  using match_reference_type = typename std::iterator_traits<IterB>::reference;
#if defined(__cpp_lib_is_invocable) && __cpp_lib_is_invocable >= 201703
  using score_type =
      std::invoke_result_t<std::decay_t<ScoreFun>, ref_reference_type,
                           match_reference_type>;
#else
  using score_type = std::result_of_t<std::decay_t<ScoreFun>(
      ref_reference_type, match_reference_type)>;
#endif

  return getBestMatchingIndices(std::move(refBegin), std::move(refEnd),
                                std::move(matchBegin), std::move(matchEnd),
                                std::forward<ScoreFun>(scoreFun),
                                std::less<score_type>());
}

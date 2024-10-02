#pragma once

#include "ringmap_matrix.hpp"
#include "ringmap_matrix_row.hpp"
#include "utils.hpp"

#include <memory>

template <typename Iterable>
void RingmapMatrix::addRead(
    Iterable &&iterable,
    std::enable_if_t<not std::is_same<std::decay_t<Iterable>,
                                      MutationMapTranscriptRead>::value>
        *) noexcept(false) {
  while (readsCount >= data.size())
    data.emplace_back();
  auto &read = data[readsCount++];
  unsigned baseIndex = 0;
  for (auto &&element : iterable) {
    if (element == 1)
      read.emplace_back(baseIndex);
    ++baseIndex;
  }
}

template <typename TranscriptRead>
void RingmapMatrix::addRead(
    TranscriptRead &&transcriptRead,
    std::enable_if_t<std::is_same<std::decay_t<TranscriptRead>,
                                  MutationMapTranscriptRead>::value>
        *) noexcept(false) {
  assert(std::ranges::is_sorted(transcriptRead.indices));
  while (readsCount >= data.size())
    data.emplace_back();
  data[readsCount] = std::forward<TranscriptRead>(transcriptRead).indices;
  data[readsCount].begin_index = transcriptRead.begin;
  data[readsCount].end_index = transcriptRead.end;
  ++readsCount;
}

template <typename Iterable>
void RingmapMatrix::keepOnlyIndices(Iterable &&iterable) noexcept(false) {
  using iterable_type = std::decay_t<Iterable>;
  const iterable_type *sortedIterable;
  std::unique_ptr<iterable_type> localIterable;
  if (std::ranges::is_sorted(iterable))
    sortedIterable = &iterable;
  else {
    localIterable.reset(new iterable_type(std::forward<Iterable>(iterable)));
    std::ranges::sort(*localIterable);
    sortedIterable = localIterable.get();
  }

  auto end = std::ranges::next(std::ranges::begin(data), readsCount);
  for (unsigned index = readsCount;;) {
    if (index == 0)
      break;
    --index;

    if (not std::ranges::binary_search(*sortedIterable, index))
      data[index] = std::move(*--end);
  }
  readsCount = static_cast<unsigned>(
      std::ranges::distance(std::ranges::begin(data), end));
  data.resize(readsCount);
}

template <typename Weights>
arma::mat RingmapMatrix::covariance(Weights &&baseWeights) const
    noexcept(false) {
  RingmapMatrix transposed = t();
  assert(transposed.t() == *this);
  arma::mat out(bases, bases, arma::fill::zeros);
  for (unsigned row = 0; row < bases; ++row) {
    double rowBaseWeight = baseWeights[row];
    if (rowBaseWeight == 0)
      continue;

    const auto &rowVector = transposed.data[row];
    for (unsigned col = row; col < bases; ++col) {
      double colBaseWeight = baseWeights[col];
      if (colBaseWeight == 0)
        continue;

      const auto &colVector = transposed.data[col];
      std::size_t count = count_intersections(
          std::ranges::begin(rowVector), std::ranges::end(rowVector),
          std::ranges::begin(colVector), std::ranges::end(colVector));

#ifndef NDEBUG
      if (row == col)
        assert(count == colVector.size());
#endif
      out(row, col) =
          static_cast<double>(count) / std::sqrt(rowBaseWeight * colBaseWeight);
    }
  }

  return arma::symmatu(out);
}

static_assert(
    std::random_access_iterator<typename RingmapMatrix::row_iterator>);
static_assert(
    std::random_access_iterator<typename RingmapMatrix::const_row_iterator>);
static_assert(
    std::random_access_iterator<typename RingmapMatrix::col_iterator>);
static_assert(
    std::random_access_iterator<typename RingmapMatrix::const_col_iterator>);
static_assert(std::random_access_iterator<
              typename RingmapMatrix::row_iterator_helper::iterator>);
static_assert(std::random_access_iterator<
              typename RingmapMatrix::const_row_iterator_helper::iterator>);
static_assert(std::random_access_iterator<
              typename RingmapMatrix::col_iterator_helper::iterator>);
static_assert(std::random_access_iterator<
              typename RingmapMatrix::const_col_iterator_helper::iterator>);

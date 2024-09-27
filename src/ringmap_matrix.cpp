#include "ringmap_matrix.hpp"
#include "ringmap_matrix_col_iterator.hpp"
#include "ringmap_matrix_row_iterator.hpp"
#include "utils.hpp"

#include <cassert>
#include <range/v3/algorithm.hpp>
#include <range/v3/numeric.hpp>

RingmapMatrix::RingmapMatrix(unsigned nBases) noexcept(false)
    : bases(nBases), readsCount(0) {}

RingmapMatrix::RingmapMatrix(unsigned nReads, unsigned nBases) noexcept(false)
    : bases(nBases), readsCount(0), data(nReads) {}

void RingmapMatrix::addModifiedIndicesRow(row_type const &row) noexcept(false) {
  assert(row.end_index >= row.begin_index);
  assert(row.end_index - row.begin_index <= bases);
  data.emplace_back(row);
  ++readsCount;
}

void RingmapMatrix::addModifiedIndicesRow(row_type &&row) noexcept(false) {
  assert(row.end_index >= row.begin_index);
  assert(row.end_index - row.begin_index <= bases);
  data.emplace_back(std::move(row));
  ++readsCount;
}

arma::vec RingmapMatrix::mean(unsigned char axis) const noexcept(false) {
  unsigned divider;
  if (axis == 0)
    divider = readsCount;
  else
    divider = bases;

  return sum(axis) / divider;
}

arma::vec RingmapMatrix::sum(unsigned char axis) const noexcept(false) {
  if (axis == 0) {
    arma::vec out(bases, arma::fill::zeros);
    for (const auto &readData : data) {
      for (unsigned index : readData)
        out[index] += 1.;
    }

    return out;
  } else {
    arma::vec out(readsCount);
    unsigned readIndex = 0;
    for (const auto &readData : data)
      out[readIndex++] = static_cast<double>(readData.size());

    return out;
  }
}

auto RingmapMatrix::operator()(unsigned row,
                               unsigned col) noexcept -> accessor_type {
  readsCount = std::max(row + 1, readsCount);
  return {data.data() + row, col};
}

auto RingmapMatrix::operator()(unsigned row, unsigned col) const noexcept
    -> const_accessor_type {
  return {data.data() + row, col};
}

void RingmapMatrix::remove_rows(unsigned begin, unsigned end) noexcept(false) {
  readsCount -= std::min(readsCount, end) - std::min(readsCount, begin);
  data.erase(ranges::next(ranges::begin(data), begin),
             ranges::next(ranges::begin(data), end));
}

void RingmapMatrix::remove_cols(unsigned begin, unsigned end) noexcept {
  bases -= std::min(bases, end) - std::min(bases, begin);
  for (auto &row : data) {
    row.erase(ranges::remove_if(row,
                                [begin, end](unsigned index) {
                                  return index >= begin and index < end;
                                }),
              ranges::end(row));
  }
}

unsigned RingmapMatrix::rows_size() const noexcept {
  assert(data.size() <= std::numeric_limits<unsigned>::max());
  return static_cast<unsigned>(data.size());
}

unsigned RingmapMatrix::cols_size() const noexcept { return bases; }

unsigned RingmapMatrix::storedReads() const noexcept { return readsCount; }

auto RingmapMatrix::row(unsigned index) noexcept -> row_accessor {
  readsCount = std::max(index + 1, readsCount);
  return {*this, data.data() + index};
}

auto RingmapMatrix::row(unsigned index) const noexcept -> const_row_accessor {
  return {*this, data.data() + index};
}

auto RingmapMatrix::col(unsigned index) noexcept -> col_accessor {
  return {*this, index};
}

auto RingmapMatrix::col(unsigned index) const noexcept -> const_col_accessor {
  return {*this, index};
}

void RingmapMatrix::shrink() noexcept { data.resize(readsCount); }

auto RingmapMatrix::rows() noexcept -> row_iterator_helper { return {*this}; }

auto RingmapMatrix::rows() const noexcept -> const_row_iterator_helper {
  return {*this};
}
auto RingmapMatrix::cols() noexcept -> col_iterator_helper { return {*this}; }

auto RingmapMatrix::cols() const noexcept -> const_col_iterator_helper {
  return {*this};
}

RingmapMatrix RingmapMatrix::t() const noexcept(false) {
  RingmapMatrix out(bases, readsCount);
  unsigned readIndex = 0;
  for (const auto &read : data) {
    for (unsigned index : read)
      out.data[index].emplace_back(readIndex);
    ++readIndex;
  }
  out.readsCount = bases;

  return out;
}

arma::mat RingmapMatrix::covariance() const noexcept(false) {
  RingmapMatrix transposed = t();
  assert(transposed.t() == *this);
  arma::mat out(bases, bases, arma::fill::zeros);
  for (unsigned row = 0; row < bases; ++row) {
    const auto &rowVector = transposed.data[row];
    for (unsigned col = row; col < bases; ++col) {
      const auto &colVector = transposed.data[col];
      unsigned count = static_cast<unsigned>(count_intersections(
          ranges::begin(rowVector), ranges::end(rowVector),
          ranges::begin(colVector), ranges::end(colVector)));

#ifndef NDEBUG
      if (row == col)
        assert(count == colVector.size());
#endif

      out(row, col) = count;
    }
  }

  return arma::symmatu(out);
}

bool RingmapMatrix::operator==(const RingmapMatrix &other) const noexcept {
  if (bases != other.bases or readsCount != other.readsCount)
    return false;

  return ranges::equal(data | ranges::view::take(readsCount),
                       other.data | ranges::view::take(readsCount));
}

bool RingmapMatrix::operator!=(const RingmapMatrix &other) const noexcept {
  if (bases != other.bases or readsCount != other.readsCount)
    return true;

  auto const end = ranges::next(ranges::begin(data), readsCount);
  auto mismatch = ranges::mismatch(
      ranges::begin(data), end, ranges::begin(other.data),
      ranges::next(ranges::begin(other.data), other.readsCount));
  return mismatch.in1 != end;
}

void RingmapMatrix::append(const RingmapMatrix &other) noexcept(false) {
  assert(bases == other.bases);
  data.resize(readsCount + other.readsCount);
  ranges::copy(other.data | ranges::view::take(other.readsCount),
               ranges::next(ranges::begin(data), readsCount));
  readsCount += other.readsCount;
}

void RingmapMatrix::shuffle() noexcept(std::is_nothrow_swappable_v<row_type>) {
  std::random_device randomDevice;
  std::mt19937 randomGenerator(randomDevice());
  ranges::shuffle(data, randomGenerator);
}

void RingmapMatrix::resize(unsigned size) noexcept(false) {
  readsCount = size;
  data.resize(readsCount);
}

auto RingmapMatrix::getIndices(unsigned rowIndex) const noexcept
    -> const row_type & {
  return data[rowIndex];
}

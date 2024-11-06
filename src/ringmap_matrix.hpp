#pragma once

#include "mutation_map_transcript_read.hpp"
#include "ringmap_matrix_accessor.hpp"
#include "ringmap_matrix_col_accessor.hpp"
#include "ringmap_matrix_col_iterator.hpp"
#include "ringmap_matrix_iterator.hpp"
#include "ringmap_matrix_iterator_helper.hpp"
#include "ringmap_matrix_row_accessor.hpp"
#include "ringmap_matrix_row_iterator.hpp"
#include "ringmap_matrix_traits.hpp"

#include <vector>

#include <armadillo>

class RingmapMatrix {
  template <typename> friend class RingmapMatrixRowIterator;
  template <typename> friend class RingmapMatrixRowAccessor;
  template <typename> friend class RingmapMatrixColIterator;
  template <typename> friend class RingmapMatrixColAccessor;
  template <typename, RingmapMatrixIteratorType>
  friend class RingmapMatrixIterator;
  template <typename> friend class RingmapMatrixAccessor;
  template <typename> friend class detail::RingmapMatrixRowIteratorHelper;
  template <typename> friend class detail::RingmapMatrixColIteratorHelper;

public:
  using value_type = ringmap_matrix::value_type;
  using row_type = ringmap_matrix::row_type;
  using matrix_type = ringmap_matrix::matrix_type;
  using accessor_type = RingmapMatrixAccessor<RingmapMatrix>;
  using const_accessor_type = RingmapMatrixAccessor<const RingmapMatrix>;
  using row_iterator = RingmapMatrixRowIterator<RingmapMatrix>;
  using const_row_iterator = RingmapMatrixRowIterator<const RingmapMatrix>;
  using row_accessor = RingmapMatrixRowAccessor<RingmapMatrix>;
  using const_row_accessor = RingmapMatrixRowAccessor<const RingmapMatrix>;
  using col_iterator = RingmapMatrixColIterator<RingmapMatrix>;
  using const_col_iterator = RingmapMatrixColIterator<const RingmapMatrix>;
  using col_accessor = RingmapMatrixColAccessor<RingmapMatrix>;
  using const_col_accessor = RingmapMatrixColAccessor<const RingmapMatrix>;
  using row_iterator_helper =
      detail::RingmapMatrixRowIteratorHelper<RingmapMatrix>;
  using const_row_iterator_helper =
      detail::RingmapMatrixRowIteratorHelper<const RingmapMatrix>;
  using col_iterator_helper =
      detail::RingmapMatrixColIteratorHelper<RingmapMatrix>;
  using const_col_iterator_helper =
      detail::RingmapMatrixColIteratorHelper<const RingmapMatrix>;

  RingmapMatrix() = default;
  explicit RingmapMatrix(unsigned nBases) noexcept(false);
  RingmapMatrix(unsigned nReads, unsigned nBases) noexcept(false);

  accessor_type operator()(unsigned row, unsigned col) noexcept;
  const_accessor_type operator()(unsigned row, unsigned col) const noexcept;

  template <typename Iterable>
  void addRead(
      Iterable &&iterable,
      std::enable_if_t<not std::is_same<std::decay_t<Iterable>,
                                        MutationMapTranscriptRead>::value> * =
          nullptr) noexcept(false);
  template <typename TranscriptRead>
  void
  addRead(TranscriptRead &&transcriptRead,
          std::enable_if_t<std::is_same<std::decay_t<TranscriptRead>,
                                        MutationMapTranscriptRead>::value> * =
              nullptr) noexcept(false);

  void addModifiedIndicesRow(row_type const &row) noexcept(false);
  void addModifiedIndicesRow(row_type &&row) noexcept(false);

  unsigned rows_size() const noexcept;
  unsigned cols_size() const noexcept;
  unsigned storedReads() const noexcept;

  row_accessor row(unsigned index) noexcept;
  const_row_accessor row(unsigned index) const noexcept;
  col_accessor col(unsigned index) noexcept;
  const_col_accessor col(unsigned index) const noexcept;

  row_iterator_helper rows() noexcept;
  const_row_iterator_helper rows() const noexcept;
  col_iterator_helper cols() noexcept;
  const_col_iterator_helper cols() const noexcept;

  bool operator==(const RingmapMatrix &other) const noexcept;
  bool operator!=(const RingmapMatrix &other) const noexcept;
  arma::vec mean(unsigned char axis = 0) const noexcept(false);
  arma::vec sum(unsigned char axis = 0) const noexcept(false);
  RingmapMatrix t() const noexcept(false);
  arma::mat covariance() const noexcept(false);
  template <typename Weights>
  arma::mat covariance(Weights &&baseWeights) const noexcept(false);

  void remove_rows(unsigned begin, unsigned end) noexcept(false);
  void remove_cols(unsigned begin, unsigned end) noexcept;
  void shrink() noexcept;
  void shuffle() noexcept(std::is_nothrow_swappable_v<row_type>);
  void resize(unsigned size) noexcept(false);
  const row_type &getIndices(unsigned rowIndex) const noexcept;

  template <typename Iterable>
  void keepOnlyIndices(Iterable &&iterable) noexcept(false);

  void append(const RingmapMatrix &other) noexcept(false);

  bool has_same_indices(RingmapMatrix const &other) const noexcept;

private:
  unsigned bases;
  unsigned readsCount;
  matrix_type data;
};

#include "ringmap_matrix_impl.hpp"

#include "ringmap_data.hpp"
#include "ringmap_matrix_row.hpp"

#include <array>
#include <cassert>
#include <random>

constexpr std::size_t n_reads = 50000;
constexpr std::size_t read_size = 75;
constexpr std::size_t window_size = 67;
constexpr std::size_t sequence_length = 200;
constexpr std::array<double, 4> modification_cum_probabilities{0.5, 0.75, 0.90,
                                                               0.97};
constexpr std::array<char, 4> valid_bases{'A', 'C', 'T', 'G'};

static std::tuple<std::vector<RingmapMatrixRow>, RingmapData>
generate_random_ringmap() {
  thread_local std::mt19937 random_gen(std::random_device{}());
  thread_local std::uniform_int_distribution<unsigned> read_begin_dist(
      std::size_t(0), sequence_length - read_size - 1);
  thread_local std::uniform_real_distribution<double> norm_dist;
  thread_local std::uniform_int_distribution<std::size_t> valid_base_index_dist(
      std::size_t(0), std::size_t(3));

  std::string sequence(sequence_length, 'N');
  std::generate(std::begin(sequence), std::end(sequence),
                [&] { return valid_bases[valid_base_index_dist(random_gen)]; });

  std::tuple<std::vector<RingmapMatrixRow>, RingmapData> out{
      std::vector<RingmapMatrixRow>(n_reads, RingmapMatrixRow(read_size)), {}};

  auto &rows = std::get<0>(out);
  auto &ringmap_data = std::get<1>(out);
  test::RingmapData::m_data(ringmap_data) = RingmapMatrix(sequence_length);
  test::RingmapData::sequence(ringmap_data) = sequence;
  test::RingmapData::startIndex(ringmap_data) = 0;
  test::RingmapData::endIndex(ringmap_data) = sequence_length;
  test::RingmapData::nSetReads(ringmap_data) = sequence_length;

  std::vector<std::size_t> targetable_bases;
  targetable_bases.reserve(sequence_length);
  {
    std::size_t base_index = 0;
    for (char base : sequence) {
      if (base == 'A' or base == 'C')
        targetable_bases.emplace_back(base_index);
      ++base_index;
    }
  }

  for (auto &&row : rows) {
    row.set_begin_index(read_begin_dist(random_gen));
    row.set_end_index(row.begin_index() + read_size);

    auto const n_modifications = static_cast<std::uint8_t>(
        std::distance(std::begin(modification_cum_probabilities),
                      std::find_if(std::begin(modification_cum_probabilities),
                                   std::end(modification_cum_probabilities),
                                   [probability = norm_dist(random_gen)](
                                       double cur_probability) {
                                     return probability <= cur_probability;
                                   })));

    row.resize(n_modifications);
    std::sample(std::begin(targetable_bases), std::end(targetable_bases),
                std::begin(row), n_modifications, random_gen);
    std::sort(std::begin(row), std::end(row));
    test::RingmapData::m_data(ringmap_data).addModifiedIndicesRow(row);
  }

  return out;
}

static std::tuple<std::vector<RingmapMatrixRow>, RingmapData>
get_window(RingmapData const &ringmap_data,
           std::vector<RingmapMatrixRow> const &rows, unsigned start_base) {
  unsigned end_base = start_base + window_size;
  std::tuple<std::vector<RingmapMatrixRow>, RingmapData> out{
      std::vector<RingmapMatrixRow>{},
      ringmap_data.get_new_range(start_base, end_base)};

  auto &filtered_rows = std::get<0>(out);

  filtered_rows.reserve(n_reads);
  for (auto &&row : rows) {
    if (row.begin_index() <= start_base and row.end_index() >= end_base) {
      RingmapMatrixRow new_row;
      new_row.set_begin_index(start_base);
      new_row.set_end_index(end_base);

      new_row.reserve(row.size());
      for (auto modified_index : row) {
        if (modified_index >= start_base and modified_index < end_base)
          new_row.emplace_back(modified_index - start_base);
      }
      filtered_rows.emplace_back(std::move(new_row));
    }
  }

  return out;
}

static void test_get_window() {
  constexpr std::size_t n_tests = 10;
  for (std::size_t test_index = 0; test_index < n_tests; ++test_index) {
    auto &&[rows, ringmap_data] = generate_random_ringmap();

    for (unsigned start_base = 0; start_base < sequence_length - window_size;
         ++start_base) {

      auto &&[filtered_rows, ringmap_window] =
          get_window(ringmap_data, rows, start_base);

      assert(test::RingmapData::m_data(ringmap_window).rows_size() ==
             filtered_rows.size());
      for (unsigned row_index = 0; row_index < filtered_rows.size();
           ++row_index) {
        assert(
            test::RingmapData::m_data(ringmap_window).getIndices(row_index) ==
            filtered_rows[row_index]);
      }
    }
  }
}

static void test_filter_window() {
  constexpr std::size_t n_tests = 10;
  for (std::size_t test_index = 0; test_index < n_tests; ++test_index) {
    auto &&[rows, ringmap_data] = generate_random_ringmap();

    unsigned const minimumModificationsPerRead =
        test::RingmapData::minimumModificationsPerRead(ringmap_data);
    const auto minimumCoverage =
        test::RingmapData::minimumCoverage(ringmap_data);
    const auto minimumModificationsPerBaseFraction =
        test::RingmapData::minimumModificationsPerBaseFraction(ringmap_data);
    std::string const sequence = ringmap_data.getSequence();

    for (unsigned start_base = 0; start_base < sequence_length - window_size;
         start_base += 10) {

      auto &&[filtered_rows, ringmap_window] =
          get_window(ringmap_data, rows, start_base);

      {
        auto const minimumModificationsPerBase =
            test::RingmapData::minimumModificationsPerBase(ringmap_window);

        std::vector<unsigned> bases_to_keep;
        for (unsigned base_index = 0; base_index < window_size; ++base_index) {
          if (sequence[base_index + start_base] != 'A' and
              sequence[base_index + start_base] != 'C')
            continue;

          const auto baseCoverage =
              test::RingmapData::baseCoverages(ringmap_window)[base_index];
          if (baseCoverage < minimumCoverage) {
            continue;
          }

          if (filtered_rows.size() <
              test::RingmapData::minimumCoverage(ringmap_window))
            continue;

          auto const modificationsOnCol = static_cast<double>(
              test::RingmapData::m_data(ringmap_window).col(base_index).sum());
          if (modificationsOnCol / baseCoverage <=
              minimumModificationsPerBaseFraction) {
            continue;
          }

          auto const base_modifications = static_cast<std::size_t>(
              std::count_if(std::begin(filtered_rows), std::end(filtered_rows),
                            [&](auto &&row) {
                              return std::binary_search(
                                  std::begin(row), std::end(row), base_index);
                            }));

          if (base_modifications <= minimumModificationsPerBase)
            continue;

          if (static_cast<double>(base_modifications) /
                  static_cast<double>(filtered_rows.size()) >
              minimumModificationsPerBaseFraction)
            bases_to_keep.emplace_back(base_index);
        }

        ringmap_window.filterBases();

        {
          auto const &old_cols_to_new =
              ringmap_window.getNonFilteredToFilteredMap();
          std::vector<unsigned> old_cols(old_cols_to_new.size());
          std::transform(std::begin(old_cols_to_new), std::end(old_cols_to_new),
                         std::begin(old_cols), [](auto &&old_col_to_new) {
                           return std::get<0>(old_col_to_new);
                         });
          std::sort(std::begin(old_cols), std::end(old_cols));
          assert(std::equal(std::begin(bases_to_keep), std::end(bases_to_keep),
                            std::begin(old_cols), std::end(old_cols)));
        }

        ringmap_window.filterReads();

        std::vector<RingmapMatrixRow> new_filtered_rows;
        new_filtered_rows.reserve(filtered_rows.size());

        for (auto &&row : filtered_rows) {
          RingmapMatrixRow new_row;
          new_row.copy_begin_end_indices(row);

          std::copy_if(std::begin(row), std::end(row),
                       std::back_inserter(new_row), [&](auto base_index) {
                         return std::binary_search(std::begin(bases_to_keep),
                                                   std::end(bases_to_keep),
                                                   base_index);
                       });
          if (new_row.size() >= minimumModificationsPerRead)
            new_filtered_rows.emplace_back(std::move(new_row));
        }

        filtered_rows = std::move(new_filtered_rows);
      }

      auto const new_cols_to_old = ringmap_window.getFilteredToNonFilteredMap();
      for (unsigned row_index = 0; row_index < filtered_rows.size();
           ++row_index) {
        auto const &ringmap_row =
            test::RingmapData::m_data(ringmap_window).getIndices(row_index);
        auto const &filtered_row = filtered_rows[row_index];
        assert(std::equal(std::begin(ringmap_row), std::end(ringmap_row),
                          std::begin(filtered_row), std::end(filtered_row),
                          [&](auto &&ringmap_col_index, auto &&base_index) {
                            auto col_iter =
                                new_cols_to_old.find(ringmap_col_index);
                            assert(col_iter != std::end(new_cols_to_old));
                            return col_iter->second == base_index;
                          }));
      }
    }
  }
}

static void test_window_covariance() {
  auto ringmap_data = std::get<1>(generate_random_ringmap());
  std::string const sequence = ringmap_data.getSequence();

  for (unsigned start_base = 0; start_base < sequence_length - window_size;
       start_base += 10) {

    auto ringmap_window =
        ringmap_data.get_new_range(start_base, start_base + window_size);
    ringmap_window.filterBases();
    ringmap_window.filterReads();

    arma::mat cov_mat;
    {
      auto const &ringmap_mat = ringmap_window.data();
      arma::mat data =
          arma::zeros(ringmap_mat.rows_size(), ringmap_mat.cols_size());
      for (unsigned row_index = 0; row_index < ringmap_mat.rows_size();
           ++row_index) {
        auto const &row_indices = ringmap_mat.getIndices(row_index);
        for (auto col_index : row_indices)
          data(row_index, col_index) = 1.;
      }

      cov_mat = data.t() * data;
    }

    assert(arma::approx_equal(cov_mat, ringmap_window.data().covariance(),
                              "absdiff", 0.001));
  }
}

int main() {
  test_get_window();
  test_filter_window();
  test_window_covariance();
}

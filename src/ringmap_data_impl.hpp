#pragma once

#include "ringmap_data.hpp"
#include <algorithm>
#include <cassert>
#include <iterator>
#include <ranges>
#include <stdexcept>

template <typename Iter>
RingmapData::RingmapData(std::string_view sequence, unsigned nReads,
                         Iter readsBegin, Iter readsEnd, Args const &args)
    : startIndex(0), endIndex(static_cast<unsigned>(sequence.size())),
      nSetReads(nReads), minimumCoverage(args.minimum_base_coverage()),
      minimumModificationsPerBase(args.minimum_modifications_per_base()),
      minimumModificationsPerRead(args.minimum_modifications_per_read()),
      minimumModificationsPerBaseFraction(
          args.minimum_modifications_per_base_fraction()),
      sequence(sequence), baseCoverages(endIndex, 0), m_data(nReads, endIndex),
      shape(args.shape()) {
  addReads(std::move(readsBegin), std::move(readsEnd));
}

template <typename Iter> void RingmapData::addReads(Iter begin, Iter end) {
#ifndef NDEBUG
  std::size_t readsCount = 0;
#endif
  std::ranges::for_each(begin, end, [&, this](auto &&read) {
    auto beginCurrentBaseCoverages =
        std::ranges::next(std::ranges::begin(baseCoverages), read.begin);
    std::ranges::transform(
        beginCurrentBaseCoverages,
        std::ranges::next(std::ranges::begin(baseCoverages), read.end),
        beginCurrentBaseCoverages, [](unsigned value) { return value + 1; });
    m_data.addRead(read);

#ifndef NDEBUG
    ++readsCount;
#endif
  });

  assert(readsCount == nSetReads);
}

template <typename Iterable>
void RingmapData::keepOnlyIndices(Iterable &&iterable) {
  m_data.keepOnlyIndices(std::forward<Iterable>(iterable));
  nSetReads = m_data.storedReads();
}

template <typename T>
std::enable_if_t<std::is_same<std::decay_t<T>, RingmapData>::value, RingmapData>
RingmapData::from(T &&other) {
  RingmapData ringmap{std::forward<T>(other).sequence,
                      data_type(std::forward<T>(other).m_data),
                      other.startIndex, other.endIndex};
  ringmap.baseCoverages = std::forward<T>(other).baseCoverages;
  ringmap.minimumCoverage = other.minimumCoverage;
  ringmap.minimumModificationsPerBase = other.minimumModificationsPerBase;
  ringmap.minimumModificationsPerRead = other.minimumModificationsPerRead;
  ringmap.minimumModificationsPerBaseFraction =
      other.minimumModificationsPerBaseFraction;

#ifndef NDEBUG
  ringmap.allowedMismatches = other.allowedMismatches;
#endif
  ringmap.shape = other.shape;

  return ringmap;
}

template <typename R>
  requires std::ranges::range<R> and
           std::same_as<std::ranges::range_reference_t<R>, RingmapData &> and
           std::random_access_iterator<std::ranges::iterator_t<R>>
void RingmapData::filter_bases_on_replicates(R &ringmap_data_range) {
  filter_bases_on_replicates_inner(
      ringmap_data_range,
      [&](auto &first_ringmap_data, auto original_base_index) {
        auto const minimum_modifications_per_base_fraction =
            first_ringmap_data.minimumModificationsPerBaseFraction;
        auto const minimum_modifications_per_base =
            first_ringmap_data.minimumModificationsPerBase;
        auto const bases_filtered = first_ringmap_data.basesFiltered;
        auto const minimum_coverage = first_ringmap_data.minimumCoverage;

        return std::ranges::all_of(
            ringmap_data_range, [&](const auto &ringmap_data) {
              auto const modifications_on_col = static_cast<double>(
                  ringmap_data.m_data
                      .col(static_cast<unsigned>(original_base_index))
                      .sum());

              return modifications_on_col >= minimum_modifications_per_base and
                     (bases_filtered or
                      (ringmap_data.baseCoverages[original_base_index] >=
                           minimum_coverage and
                       modifications_on_col /
                               ringmap_data.baseCoverages[original_base_index] >
                           minimum_modifications_per_base_fraction));
            });
      });
}

template <typename R>
  requires std::ranges::range<R> and
           std::same_as<std::ranges::range_reference_t<R>, RingmapData &> and
           std::random_access_iterator<std::ranges::iterator_t<R>>
void RingmapData::filter_bases_on_replicates_for_assignments(
    R &ringmap_data_range) {
  filter_bases_on_replicates_inner(ringmap_data_range,
                                   [](auto const &, auto) { return true; });
}

template <typename R, typename F>
  requires std::ranges::range<R> and
           std::same_as<std::ranges::range_reference_t<R>, RingmapData &> and
           std::random_access_iterator<std::ranges::iterator_t<R>>
void RingmapData::filter_bases_on_replicates_inner(R &ringmap_data_range,
                                                   F base_filter) {
  assert(std::ranges::begin(ringmap_data_range) !=
         std::ranges::end(ringmap_data_range));
  auto const &first_ringmap_data = *std::ranges::begin(ringmap_data_range);
  assert(std::ranges::all_of(
      ringmap_data_range | std::views::drop(1), [&](auto const &ringmap_data) {
        return ringmap_data.m_data.cols_size() ==
                   first_ringmap_data.m_data.cols_size() and
               ringmap_data.basesMask == first_ringmap_data.basesMask and
               ringmap_data.basesFiltered ==
                   first_ringmap_data.basesFiltered and
               ringmap_data.sequence == first_ringmap_data.sequence and
               ringmap_data.shape == first_ringmap_data.shape and
               ringmap_data.minimumModificationsPerBase ==
                   first_ringmap_data.minimumModificationsPerBase and
               ringmap_data.minimumCoverage ==
                   first_ringmap_data.minimumCoverage and
               ringmap_data.minimumModificationsPerBaseFraction ==
                   first_ringmap_data.minimumModificationsPerBaseFraction and
               ringmap_data.oldColsToNew == first_ringmap_data.oldColsToNew;
      }));

  assert(std::ranges::all_of(ringmap_data_range, [&](auto const &ringmap_data) {
    return ringmap_data.basesFiltered or
           std::ranges::all_of(
               ringmap_data.baseCoverages,
               [nReads = ringmap_data.m_data.rows_size()](unsigned coverage) {
                 return coverage <= nReads;
               });
  }));

  auto allowed_bases = ([&]() {
    if (first_ringmap_data.basesMask.empty() or
        first_ringmap_data.basesFiltered) {
      return std::views::iota(0u, static_cast<unsigned>(
                                      first_ringmap_data.m_data.cols_size())) |
             std::ranges::to<std::vector>();
    } else {
      return std::views::zip(std::views::iota(0u),
                             first_ringmap_data.basesMask) |
             std::views::filter(
                 [](auto const &pair) { return std::get<1>(pair); }) |
             std::views::transform(
                 [](auto &&pair) { return std::get<0>(pair); }) |
             std::ranges::to<std::vector>();
    }
  })();

  auto filtered_bases =
      allowed_bases | std::views::transform([&](auto const base_index) {
        if (first_ringmap_data.basesFiltered) {
          auto const old_col_to_new_iter = std::ranges::find_if(
              first_ringmap_data.oldColsToNew, [&](const auto old_and_new) {
                return old_and_new.second == base_index;
              });

          assert(old_col_to_new_iter !=
                 std::ranges::end(first_ringmap_data.oldColsToNew));
          return std::make_pair(base_index, old_col_to_new_iter->first);
        } else {
          return std::make_pair(base_index, base_index);
        }
      }) |
      std::views::filter([&](auto const &pair) {
        auto const original_base_index = std::get<0>(pair);
        auto const base_index = std::get<1>(pair);
        auto const base = static_cast<char>(std::toupper(
            static_cast<int>(first_ringmap_data.sequence[base_index])));

        return (first_ringmap_data.shape or base == 'C' or base == 'A') and
               base_filter(first_ringmap_data, original_base_index);
      }) |
      std::ranges::to<std::vector>();

  if (first_ringmap_data.oldColsToNew.empty() or
      std::size(filtered_bases) != std::size(allowed_bases)) {
    std::map<unsigned, unsigned> old_cols_to_new;
    for (auto &&[new_col_index, col_index_base_index] :
         std::views::zip(std::views::iota(0u), filtered_bases)) {
      old_cols_to_new.emplace(
          static_cast<unsigned>(std::get<1>(col_index_base_index)),
          new_col_index);
    }

    auto const cols_to_use =
        filtered_bases |
        std::views::transform([](auto pair) { return std::get<0>(pair); }) |
        std::ranges::to<std::vector>();

    auto update_ringmap_data = [&](auto &ringmap_data) {
      RingmapMatrix filtered_data(
          ringmap_data.m_data.rows_size(),
          static_cast<unsigned>(std::size(cols_to_use)));
      for (auto const [filtered_col_index, col_index] :
           std::views::zip(std::views::iota(0u), cols_to_use)) {
        filtered_data.col(filtered_col_index) =
            ringmap_data.m_data.col(col_index);
      }

      for (auto &&[filtered_row, data_row] :
           std::views::zip(filtered_data.rows(), ringmap_data.m_data.rows())) {
        filtered_row.copy_begin_end_indices(data_row);
        filtered_row.copy_window_begin_end_indices(data_row);
      }

      ringmap_data.m_data = std::move(filtered_data);
      ringmap_data.baseCoverages.erase(
          std::ranges::transform(filtered_bases |
                                     std::views::transform([](auto pair) {
                                       return std::get<0>(pair);
                                     }),
                                 std::ranges::begin(ringmap_data.baseCoverages),
                                 [&](std::size_t index) {
                                   return ringmap_data.baseCoverages[index];
                                 })
              .out,
          std::ranges::end(ringmap_data.baseCoverages));
      ringmap_data.basesFiltered = true;
    };

    auto const ringmap_data_size = std::ranges::size(ringmap_data_range);
    for (auto &ringmap_data :
         ringmap_data_range | std::views::take(ringmap_data_size - 1)) {
      ringmap_data.oldColsToNew = old_cols_to_new;
      update_ringmap_data(ringmap_data);
    }

    auto &last_ringmap_data =
        *std::ranges::prev(std::ranges::end(ringmap_data_range));
    last_ringmap_data.oldColsToNew = std::move(old_cols_to_new);
    update_ringmap_data(last_ringmap_data);
  }
}

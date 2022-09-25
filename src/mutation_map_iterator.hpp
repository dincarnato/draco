#pragma once

#include <ios>
#include <type_traits>

#include <range/v3/core.hpp>

class MutationMapTranscript;

template <typename MMap> class MutationMapIterator {
  friend class MutationMap;

public:
  using value_type = MutationMapTranscript;
  using reference = std::conditional_t<std::is_const_v<MMap>,
                                       const value_type &, value_type &>;
  using pointer = std::conditional_t<std::is_const_v<MMap>, const value_type *,
                                     value_type *>;
  using iterator_category = ranges::input_iterator_tag;
  using difference_type = std::ptrdiff_t;

  MutationMapIterator() = default;
  MutationMapIterator(MMap &mutationMap, std::streampos offset = 0) noexcept;

  using transcripts_type = typename MMap::transcripts_type;
  using self = MutationMapIterator;

  reference operator*() const noexcept(false);
  pointer operator->() const noexcept(false);
  self &operator++() noexcept;
  self operator++(int) noexcept;
  difference_type operator-(const self &other) const noexcept;
  bool operator==(const self &other) const noexcept(false);
  bool operator!=(const self &other) const noexcept(false);

private:
  template <typename T = MMap>
  std::enable_if_t<not std::is_const<T>::value, bool>
  loadOtherTranscripts() const noexcept(false);

  template <typename T = MMap>
  void checkOffset(std::enable_if_t<std::is_const<T>::value> * = nullptr) const
      noexcept(false);
  template <typename T = MMap>
  void
  checkOffset(std::enable_if_t<not std::is_const<T>::value> * = nullptr) const
      noexcept(false);

  MMap *mutationMap = nullptr;
  std::streampos offset;
};

#include "mutation_map_iterator_impl.hpp"

#pragma once

#include <ranges>
#include <vector>

namespace more_ranges {

/**
 * \brief Creates a to_vector_closure for operator()
 */
struct to_vector_adapter {
  struct closure {
    /**
     * \brief Gets a vector of a given range.
     * \tparam R type of range that gets converted to a vector.
     * \param r range that gets converted to a vector.
     * \return vector from the given range.
     */
    template <std::ranges::range R> constexpr auto operator()(R &&r) const {
      std::vector<std::ranges::range_value_t<R>> v;

      // if we can get a size, reserve that much
      if constexpr (requires { std::ranges::size(r); }) {
        v.reserve(std::ranges::size(r));
      }

      v.insert(v.begin(), r.begin(), r.end());

      return v;
    }
  };

  /**
   * \brief Gets a closure to convert the range to a vector.
   * \return A to_vector_closure that will convert the range to a vector.
   */
  constexpr auto operator()() const -> closure { return closure{}; }

  template <std::ranges::range R> constexpr auto operator()(R &&r) {
    return closure{}(std::forward<R>(r));
  }
};

inline to_vector_adapter to_vector;

/**
 * \brief A range pipe that results in a vector.
 * \tparam R type of range that gets converted to a vector.
 * \param r range that gets converted to a vector.
 * \param a used to create the vector.
 * \return a vector from the given range.
 */
template <std::ranges::range R>
constexpr auto operator|(R &&r, to_vector_adapter::closure const &a) {
  return a(std::forward<R>(r));
}

} /* namespace more_ranges */

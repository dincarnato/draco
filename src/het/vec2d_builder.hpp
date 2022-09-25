#pragma once

#include "detail/vec2d_element_builder.hpp"

#include <tuple>

namespace het {

struct vec2d_builder {
  static constexpr auto build() noexcept {
    return detail::vec2d_element_builder_ct<
        void, void, std::integral_constant<std::size_t, 0>>{};
  }
};

} // namespace het

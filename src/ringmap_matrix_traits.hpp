#pragma once

#include <vector>

struct RingmapMatrixRow;

namespace ringmap_matrix {

using value_type = bool;
using base_index_type = unsigned;
using row_type = RingmapMatrixRow;
using matrix_type = std::vector<row_type>;

} // namespace ringmap_matrix

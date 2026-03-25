#include "../reassignment.hpp"

namespace test {

struct Reassignment : public ::Reassignment {
  Reassignment(::Reassignment &&reassignment)
      : ::Reassignment(std::move(reassignment)) {}

  using ::Reassignment::
      reweight_and_reassign_with_expectation_maximization_iteration;
};

} // namespace test

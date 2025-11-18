#include "clusters_replicates.hpp"

using namespace clusters_replicates;

static void test_distances_size() {
  assert(distances_size(1, 5) == 0);
  assert(distances_size(2, 5) == 25);
  assert(distances_size(4, 5) == 25 * 6);
  assert(distances_size(6, 4) == 16 * 15);
}

int main() { test_distances_size(); }

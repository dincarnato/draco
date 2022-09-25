#include <ringmap_data.hpp>

#include <cassert>

static constexpr const char *sequence = "ACCACCAACA";

int main(int argc, const char *argv[]) {
  assert(argc == 2);
  const char *const filename = argv[1];

  Args args;
  RingmapData ringmapData(filename, sequence, args);

  auto const secondPartBegin =
      static_cast<unsigned>(static_cast<double>(ringmapData.size()) / 3 * 2);
  RingmapData split1 = RingmapData::from(ringmapData);
  split1.resize(secondPartBegin);

  RingmapData split2 = RingmapData::from(ringmapData);
  split2.remove(0, secondPartBegin);

  assert(split1 + split2 == ringmapData);
}

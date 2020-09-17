#include <ringmap_data.hpp>

#include <random>
#include <cassert>

static constexpr const char* sequence = "ACCACCAACA";
static constexpr unsigned loops = 1000;

int main(int argc, const char* argv[])
{
    assert(argc == 2);
    const char* const filename = argv[1];
    Args args;

    RingmapData ringmapData(filename, sequence, args);

    std::random_device randomDevice;
    std::mt19937 randomGenerator(randomDevice());

    /* Row-wise shuffle */
    {
        std::vector<unsigned> sums;
        for(auto&& row : ringmapData.data().rows())
            sums.emplace_back(row.sum());

        for(unsigned loop = 0; loop < loops; ++loop)
        {
            auto matrix = ringmapData.data();
            for(auto&& row : matrix.rows())
                row.shuffle(randomGenerator);

            assert(matrix != ringmapData.data());

            auto sumIter = std::begin(sums);
            for(const auto& row : matrix.rows())
                assert(row.sum() == *sumIter++);
        }
    }

    /* Col-wise shuffle */
    {
        std::vector<unsigned> sums;
        for(auto&& col : ringmapData.data().cols())
            sums.emplace_back(col.sum());

        for(unsigned loop = 0; loop < loops; ++loop)
        {
            auto matrix = ringmapData.data();
            for(auto&& col : matrix.cols())
                col.shuffle(randomGenerator);

            assert(matrix != ringmapData.data());

            auto sumIter = std::begin(sums);
            for(const auto& col : matrix.cols())
                assert(col.sum() == *sumIter++);
        }
    }
}

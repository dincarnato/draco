#include "combinator.hpp"

#include <algorithm>
#include <array>

int main()
{
    constexpr unsigned N = 8;
    {
        Combinator combinator(8, 3);
        auto iter = std::begin(combinator);
        for(unsigned a = 0; a < N; ++a)
        {
            for(unsigned b = a + 1; b < N; ++b)
            {
                for(unsigned c = b + 1; c < N; ++c, ++iter)
                {
                    if(iter == std::end(combinator))
                        return -1;

                    auto&& combination = *iter;
                    std::array<unsigned, 3> realCombination{a, b, c};
                    if(not std::equal(std::begin(combination), std::end(combination), std::begin(realCombination)))
                        return -1;
                }
            }
        }

        if(iter != std::end(combinator))
            return -1;
    }

    {
        std::vector<std::size_t> values{4, 7, 10, 16, 21, 22, 30, 723};
        Combinator combinator(std::begin(values), std::end(values), 3);
        auto iter = std::begin(combinator);
        for(unsigned a = 0; a < N; ++a)
        {
            for(unsigned b = a + 1; b < N; ++b)
            {
                for(unsigned c = b + 1; c < N; ++c, ++iter)
                {
                    if(iter == std::end(combinator))
                        return -1;

                    auto&& combination = *iter;
                    std::array<std::size_t, 3> realCombination{values[a], values[b], values[c]};
                    if(not std::equal(std::begin(combination), std::end(combination), std::begin(realCombination)))
                        return -1;
                }
            }
        }

        if(iter != std::end(combinator))
            return -1;
    }
}
